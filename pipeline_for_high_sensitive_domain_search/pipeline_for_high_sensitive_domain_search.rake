

# {{{ procedures
WriteBatch  = lambda do |outs, jdir, t|
	outs.each_slice(10000).with_index(1){ |ls, idx|
		open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.#{idx}", "w"){ |fjob|
			fjob.puts ls
		}
	}
end

RunBatch    = lambda do |jdir, queue, nthreads, mem, wtime, ncpus, tdir|
	Dir["#{jdir}/*"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin|
		if queue != ""
			raise("`--queue #{queue}': invalid queue") unless %w|JP1 JP4 JP10 cdb|.include?(queue)
			sh "qsubarraywww -q #{queue} -l ncpus=#{nthreads} -l mem=#{mem}gb -l walltime=#{wtime} #{fin}"
		elsif ncpus != ""
			raise("`--ncpus #{ncpus}': not an integer") if ncpus !~ /^\d+$/
			sh "parallel --tmpdir #{tdir} --jobs #{ncpus} <#{fin}" ### use tmpdir for temporary starge
		else
			sh "sh #{fin}"
		end
	}
end

PrintStatus = lambda do |current, total, status, t|
	puts ""
	puts "\e[1;32m===== #{Time.now}\e[0m"
	puts "\e[1;32m===== step #{current} / #{total} (#{t.name}) -- #{status}\e[0m"
	puts ""
	$stdout.flush
end

CheckVersion = lambda do |commands|
	commands.each{ |command|
		str = case command
					when "ruby"
						%{ruby --version 2>&1}
					when "makeblastdb"
						%{makeblastdb -version 2>&1}
					when "tblastx"
						%{tblastx -version 2>&1}
					when "R"
						%{LANG=C R --version 2>&1}
					when "ape"
						%{LANG=C R --quiet --no-save --no-restore -e "packageVersion('ape')" 2>&1}
					when "phangorn"
						%{LANG=C R --quiet --no-save --no-restore -e "packageVersion('phangorn')" 2>&1}
					when "hhsearch"
						%{hhsearch 2>&1 |grep "^HHsearch"}
					when "jackhmmer"
						%{jackhmmer -h 2>&1 |head -2}
					when "reformat.pl"
						%{reformat.pl 2>&1 |sed -ne '2p'}
					end
		puts ""
		puts "\e[1;32m===== check version: #{command}\e[0m"
		puts ""
		puts "$ #{str}"
		### run
		puts `#{str}`
		### flush
		$stdout.flush
	}
end
# }}} procedures


# {{{ default (run all tasks)
task :default do
	### check version
	CheckVersion.call(%w|jackhmmer hhsearch reformat.pl ruby|)

	### tasks
	tasks = %w|01-1.prep_query_files 01-2.parse_hhdb_ent 01-3.jackhmmer 01-4.reformat.pl 01-5.hhsearch 01-6.parse_result|

	In       = ENV["in"]
	Odir     = ENV["dir"]
	Jackdb   = ENV["jackdb"]
	Hhdb     = ENV["hhdb"]
	
	HhdbData = "#{Hhdb}_hhm.ffdata" # to parse id-acc-def
	raise("\e[1;31mError:\e[0m #{HhdbData} should exist but not found") unless File.exist?(HhdbData)

	Iter     = ENV["iter"].to_i    # default: 5
	IncE     = ENV["incE"].to_f    # default: 0.001
	IncdomE  = ENV["incdomE"].to_f # default: 0.001
	# E        = ENV["e"].to_f       # default: 0.001, for hhsearch inclusion threshold

	Nreport  = ENV["nreport"].to_i # default: 3
	Evalue   = ENV["evalue"].to_f  # default: 10
	Pvalue   = ENV["pvalue"].to_f  # default: 1
	Prob     = ENV["prob"].to_f    # default: 30

	Nthreads = "1".to_i              
	Mem      = Nthreads * 12
	Qname    = ENV["queue"]||""
	Wtime    = ENV["wtime"]||"24:00:00"
	Ncpus    = ENV["ncpus"]||""    

  Tdir     = "#{Odir}/tmp"; mkdir_p Tdir

	NumStep  = tasks.size
	tasks.each.with_index(1){ |task, idx|
		Rake::Task[task].invoke(idx)
	}
end
# }}} default (run all tasks)


# {{{ tasks 01
task "01-1.prep_query_files", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	Dir1     = "#{Odir}/work"; mkdir_p Dir1
	Dir2     = "#{Odir}/info"; mkdir_p Dir2
	Fin0     = "#{Dir2}/query-stats.txt"
	# Fin1     = "#{Dir2}/query.list"
	fins     = In.split(",")  ### input file
	stats    = []
	qlist    = []
	fnames   = {}
	Fnames   = []
	Fname2pids = Hash.new{ |h, i| h[i] = [] }
	Group2pids = Hash.new{ |h, i| h[i] = [] }

	### detect whether input is fasta or list
	flag = case fins.size
				 when 1 
					ents = IO.read(fins[0]).split(/^>/)
					flag = ents.size > 1 ? "fasta" : "list" ## if input file is multiple, regard as fasta file
				 else "fasta" ## fin is fasta if multiple file is given.
				 end

	### fetch fasta files from list
	if flag == "list"
		fins = IO.readlines(fins[0]).map{ |l| l.strip }
	end

	fins.each{ |fin| ## fin would be FASTA file or list of FASTA file
		### take file name (fname)
		fname = File.basename(fin)

		### uniq file name check
		raise("\e[1;31mError:\e[0m File name should be uniq. '#{fname}' is given multiple times.") if fnames[fname]
		fnames[fname] = 1

		group2out = Hash.new{ |h, i| h[i] = {} } # group2out[group][id] = [output fasta]
		pid2len   = {}
		ents      = IO.read(fin).split(/^>/)

		ents[1..-1].each{ |ent|
			lab, *seq = ent.split("\n")
			group     = lab[/belongTo:(\S+)/, 1] || "ungrouped"
			pid       = lab.split(/\s+/)[0]
			len       = seq.join.gsub(/[^A-Za-z]/, "").size

			### uniq pid check
			raise("\e[1;31mError:\e[0m Sequence ID is not uniq. '#{pid}' is found multiple times in #{fin}.") if pid2len[pid]
			pid2len[pid] = len

			### store sequnece
			group2out[group][pid] = [">"+pid, seq.join.gsub(/[^A-Za-z]/, "")] ### The stop codon '*' must be removed. If '*' is included, reformat.pl raise error.

			### store proteins in given order
			Fname2pids[fname] << [group, pid]
			Group2pids[group] << pid
		}

		num_group         = group2out.keys.size
		num_grouped_seq   = 0
		num_ungrouped_seq = 0

		group2out.each{ |group, pid2ent|
			### calculate stats
			if group != "ungrouped"
				num_grouped_seq += pid2ent.keys.size
			else
				num_ungrouped_seq += pid2ent.keys.size
			end

			### make split fasta
			pid2ent.each{ |pid, ent|
				odir = "#{Dir1}/#{fname}/#{group}/#{pid}"; mkdir_p odir
				open("#{odir}/query.faa", "w"){ |fout| fout.puts ent }
			}
		}

		### store stats
		stats  << [fname, num_group, num_grouped_seq, num_ungrouped_seq]*"\t"
		Fnames << fname

		### write query protein list
		odir   = "#{Dir2}/query"; mkdir_p odir
		open("#{odir}/#{fname}.list", "w"){ |fout|
			group2out.each{ |group, pid2ent|
				pid2ent.each{ |pid, ent|
					fout.puts [pid, pid2len[pid], group]*"\t"
				}
			}
		}
	}

	### write stats
	open(Fin0, "w"){ |fout|
		fout.puts %w|in_file #group #grouped_protein #ungrouped_protein|*"\t"
	 	fout.puts stats
	}
end
task "01-2.parse_hhdb_ent", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	HhdbList = "#{Dir2}/hhdb.list"
	## -a option: treat binary file as text file (since `file pfam_hhm.ffdata' is data)
	sh %|grep -a "^NAME" #{HhdbData} >#{HhdbList}|
end
task "01-3.jackhmmer", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	jdir     = "#{Odir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	outs     = []

	Fnames.each{ |fname|
		Fname2pids[fname].each{ |group, pid|
			faa  = "#{Dir1}/#{fname}/#{group}/#{pid}/query.faa"
			pref = File.dirname(faa) + "/jack" ## pref = #{dir}/jack
			cmd  = "jackhmmer --cpu 1 -N #{Iter} -o #{pref}.out --tblout #{pref}.tblout --chkhmm #{pref} --chkali #{pref} \
			--incE #{IncE} --incdomE #{IncdomE} --notextw --noali #{faa} #{Jackdb}".gsub(/\s+/, " ")
			outs << cmd
		}
	}

	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, Qname, Nthreads, Mem, Wtime, Ncpus, Tdir)
end
task "01-4.reformat.pl", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	jdir     = "#{Odir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	outs     = []

	Fnames.each{ |fname|
		Fname2pids[fname].each{ |group, pid|
			fins = Dir["#{Dir1}/#{fname}/#{group}/#{pid}/jack-*.sto"]
			next if fins.size == 0

			fin  = fins.sort_by{ |fin| fin[/-(\d+)\.sto$/, 1].to_i * -1 }[0] ## select .sto of last iteration
			fout = fin.sub(/\.sto$/, ".a2m")
			log  = fin.sub(/\.sto$/, ".sto.reformat.log")

			outs << "reformat.pl #{fin} #{fout} >#{log} 2>&1"
		}
	}

	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, Qname, Nthreads, Mem, Wtime, Ncpus, Tdir)
end
task "01-5.hhsearch", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	jdir     = "#{Odir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	outs     = []

	Fnames.each{ |fname|
		Fname2pids[fname].each{ |group, pid|
			fin = Dir["#{Dir1}/#{fname}/#{group}/#{pid}/jack-*.a2m"][0]
			next unless fin
			pref = File.dirname(fin) + "/hhsearch" ## pref = #{dir}/jack
			outs << "hhsearch -d #{Hhdb} -i #{fin} -o #{pref}.hhr >#{pref}.log 2>&1"
		}
	}

	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, Qname, Nthreads, Mem, Wtime, Ncpus, Tdir)
end
desc "01-6.parse_result"
task "01-6.parse_result", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	header   = %w|protein query_len query_num_iter query_num_seq query_num_seq_used template_acc template_name template_desc
	probability e-value p-value score cols query_pos template_pos template_len rank|

	### parse Pfam desc
	acc2info = {}
	IO.readlines(HhdbList).each{ |l|
		# NAME  PF00001.20 ; 7tm_1 ; 7 transmembrane receptor (rhodopsin family)
		acc, name, desc = l.chomp.split(" ; ")
		acc = acc[/^NAME\s+(\S+)/, 1]
		acc2info[acc] = [name, desc]
	}

	### main
	Fnames.each{ |fname|
		### parse protein group/length
		fpinfo = "#{Dir2}/query/#{fname}.list" ### [pid, len, group]
		group2pinfo = Hash.new{ |h, i| h[i] = {} }
		IO.readlines(fpinfo).each{ |l|
			pid, len, group = l.chomp.split("\t")
			group2pinfo[group][pid] = len
		}

		Dir["#{Dir1}/#{fname}/*"].sort.each{ |dgroup|
			group = File.basename(dgroup)
			pinfo = group2pinfo[group]
			odir  = "#{Odir}/result/#{fname}/#{group}"; mkdir_p odir

			### store output for each group
			pid2best = {}
			pid2topN = {}

			Group2pids[group].each{ |pid|
				### parse iter/q_num_1 
				a2m     = Dir["#{dgroup}/#{pid}/jack-*.a2m"][0]
				iter    = File.basename(a2m)[/jack-(\d+)\.a2m/, 1] ## number of iteration for build queryHMM
				q_num_1 = IO.read(a2m).split(/^>/)[1..-1].size  ## number of sequence included in queryHMM

				### store [pid, len, iter]
				len  = pinfo[pid]
				out  = []
				Nreport.times do
					out << [pid, len, iter, q_num_1]
				end

				fin  = "#{dgroup}/#{pid}/hhsearch.hhr"
				if File.exist?(fin)
					ls = IO.readlines(fin)

					### parse query_num_seq (after filters such as non-redundancy filter, ...)
					q_num_2 = ls[2][/(\d+) out of \d+/, 1]
					neff    = ls[3][/^Neff\s+(\d+)/, 1]
					# out = out.map{ |a| a += [q_num_2, neff] }

					### parse best hit
					ls.each{ |l|
            ### invalid byte sequence check (require ruby version >=2.1)
            _l = l.scrub("!")
            $stderr.puts "The line inlucde invalid byte: #{l}" if _l != l
            l = _l

						#  No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
						#   1 PF04908.12 ; SH3BGR ; SH3-bind  99.3 5.9E-16 3.6E-20   88.3   0.0   75    1-75      2-91  (98)
						#   1 PF02463.18 ; SMC_N ; RecF/RecN  99.3 2.4E-16 1.4E-20  185.7   0.0  161  355-524  1100-1270(1271)
						if l =~ /^\s+(\d+)\s+(\S+)/ # top 3 hit
							rank, acc  = $1.to_i, $2
							next if rank > Nreport

							name, desc = acc2info[acc]
							scores = l[34..-1].strip.split(/\s+/)
							prob, e_val, p_val, score, ss = scores[0..4].map(&:to_f)
							cols, q_pos, t_pos, t_len     = scores[5..8]
							if t_len 
								t_len = t_len[1..-2] ## (98) --> 98
							else
								t_pos, t_len = t_pos.split("(") ## 1100-1270(1271)
								t_len        = t_len[0..-2]     ## 1271) --> 1271
							end

							if (prob > Prob) and (e_val < Evalue) and (p_val < Pvalue) ## significant topN hit
								out[rank-1] += [q_num_2, acc, name, desc, prob, e_val, p_val, score, cols, q_pos, t_pos, t_len, rank]
							end
						end
					}
				end
				out = out.map{ |a| a.size == 4 ? a += ["-"]*13 : a }
				pid2best[pid] = out[0]*"\t"
				pid2topN[pid] = out.map.with_index{ |a, idx| ## besthit --> always include (even if "-"), others --> include if significant
					idx == 0 ? a*"\t" : (a[4] == "-" ? nil : a*"\t")
				}.compact
			}

			### make output
			odir = "#{Odir}/result/#{fname}/#{group}"
			open("#{odir}/besthit.tsv", "w"){ |fout|
				fout.puts header*"\t"
				pid2best.each{ |pid, l| fout.puts l }
			}
			open("#{odir}/top#{Nreport}hit.tsv", "w"){ |fout|
				fout.puts header*"\t"
				pid2topN.each{ |pid, out| fout.puts out }
			}
		}
	}
end
# }}} tasks 01


