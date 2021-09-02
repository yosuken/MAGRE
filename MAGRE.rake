

# {{{ procedures
WriteBatch  = lambda do |t, jobdir, outs|
	jdir = "#{jobdir}/#{t.name.split(":")[-1]}"; mkdir_p jdir unless File.directory?(jdir)
  jnum = outs.size

  if jnum > 0
    outs.each_slice(jnum).with_index(1){ |ls, idx| ## always 1 file
      open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.sh", "w"){ |fjob|
        fjob.puts ls
      }
    }
  else
    open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.sh", "w"){ |fjob| } ## clear job file
  end
end

RunBatch    = lambda do |t, jobdir, ncpu, logdir|
	jdir = "#{jobdir}/#{t.name.split(":")[-1]}"
  ldir = "#{logdir}/#{t.name.split(":")[-1]}"; mkdir_p ldir unless File.directory?(ldir)

  Dir["#{jdir}/*.sh"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin| ## always 1 or 0 file
    next if File.zero?(fin)
    sh "parallel --jobs #{ncpu} --joblog #{ldir}/parallel.log <#{fin}"
    # if ncpu > 1
    #   sh "parallel --jobs #{ncpu} --joblog #{ldir}/parallel.log <#{fin}"
    # else
    # 	# sh "bash -c #{fin}" ## --> permission denied
    # 	sh "bash #{fin}"
    # end
	}
  open("#{ldir}/exit", "w"){ |fw| fw.puts "exit at #{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}" }
end

PrintStatus = lambda do |current, total, status, t|
	puts ""
	puts "\e[1;32m===== #{Time.now}\e[0m"
	puts "\e[1;32m===== step #{current} / #{total} (#{t.name}) -- #{status}\e[0m"
	puts ""
	$stdout.flush
end

CheckVersion = lambda do |commands, cmd2path|
	commands.each{ |command|
    path = cmd2path[command]||command
    f = if path.split("/").size == 1 then `which #{path}`.strip
        elsif File.exist?(path) then path
        else `which #{path}`.strip
        end
  # commands  = %w|ruby wrapper_phage_contigs_sorter_iPlant.pl prodigal jackhmmer hhsearch reformat.pl parallel|
		str = case command
					when "ruby"
						%|ruby --version 2>&1|
          when "wrapper_phage_contigs_sorter_iPlant.pl"
						%{#{path} -h 2>&1} ### version can not be checked
					when "prodigal"
						%{#{path} -v 2>&1 |head -n 2}
					when "jackhmmer"
						%{#{path} -h 2>&1 |head -n 2}
					when "hhsearch"
						%{#{path} -h 2>&1 |head -n 2}
          when "reformat.pl"
						%{#{path} -h 2>&1 |head -n 2}
					when "blastn"
						%|blastn -version 2>&1|
					when "ssearch36"
						%{ssearch36 2>&1 |head -n 7 |tail -n 3}
					when "parallel"
						%{LANG=C #{path} --version 2>&1 |head -n 1}
					# when "barrnap"
					# 	%{#{path} -h 2>&1 |head -n 2}
					# when "tRNAscan-SE"
          #   %{head -n 7 #{f}}
					# when "minced"
					# 	%|#{path} --version 2>&1|
					# when "signalp"
					# 	%|#{path} -version 2>&1|
					# when "tmhmm"
          #   %{head -n 3 #{f}} # This is version 2.0c of tmhmm
					end
		puts ""
		puts "\e[1;32m===== check version: #{command}\e[0m"
		puts ""
		puts "path: #{f}"
		puts ""
		$stdout.flush
		# puts "$ #{str}"
		### run
    out = `#{str}`
    puts out
		$stdout.flush
    raise("\e[1;31m#{str}: command not found\e[0m.") if out =~ /command not found\./
	}
end
# }}} procedures


# {{{ task controller
task :default do
	### define tasks
  tasks = []
  tasks << "01-1a.validate_input"
  tasks << "01-1b.validate_dbs"
  tasks << "01-1c.make_ctg_list"
  tasks << "01-2a.tetranucleotide"
  tasks << "01-2b.prepare_for_pca"
  tasks << "01-2c.pca"
  tasks << "01-3a.prodigal"
  tasks << "01-3b.identify_codon_table"
  tasks << "01-4.other_tasks"
  tasks << "01-5.make_contig_table"
  tasks << "01-6.make_fasta"

  # if ENV["only_detect"] == "true"
  #   Tasks = tasks[0, 5]
  # else
  #   Tasks = tasks
  # end
  Tasks = tasks

	### constants from arguments
	Odir      = ENV["outdir"]           ## output directory
	OdirExist = ENV["outdir_exist"]     ## output directory exist? ("true" or "")
	Fins      = ENV["input"]            ## input genome sequence(s)
	Fcfg      = ENV["config"]           ## config file
  Ncpu      = ENV["ncpus"].to_i       ## num of CPUs
  Fcov      = ENV["fcov"]             ## file of coverage
  Category  = ENV["category"]         ## category of remove (combination of 'taxonomy', 'coverage', 'tetranuc', 'virsorter', 'terL', 'circular')
  MinLen    = ENV["minlen"]           ## minimum contig length (in bp) to retain

	### constants
  Errmsg      = "\e[1;31mError:\e[0m"
  Warmsg      = "\e[1;35mWarning:\e[0m"
  Catemsg     = "\e[1;35m[Cateogry]\e[0m"
  Codons      = %w|4 11|
  Ncat        = [Ncpu, 4].min        ## num of cpu for each CAT/virsorter/pipeline_for_high_sensitive_domain_search
  NcatPara    = Ncpu / Ncat          ## num of parallel for each CAT/virsorter/pipeline_for_high_sensitive_domain_search
  ChunkCheckm = 20

  $valTools = {}
  keys = %w|hhsuite_pfam_db_prefix cat_bat_nr_db cat_bat_tx_db virsorter_data|
  keys.each{ |i| $valTools[i] = i }

	### parse contiguation file
  IO.readlines(Fcfg).each{ |l|
    l = l.rstrip
    case l
    when /^#/, /^$/
    else
      tool, path = l.split(/\s+/)
      $valTools[tool] = path if $valTools[tool] ## if key is already set
    end
  }
  # p $valTools

	### check version
  commands  = %w|ruby wrapper_phage_contigs_sorter_iPlant.pl prodigal jackhmmer hhsearch reformat.pl parallel|
	# commands += %w|parallel| if Ncpu > 1
	CheckVersion.call(commands, $valTools)

  ## dir path
  Ddir      = "#{File.dirname(__FILE__)}/data"
  Sdir      = "#{File.dirname(__FILE__)}/script"
  Jobdir    = "#{Odir}/batch"
  Tdir      = "#{Odir}/tmp"
  Idir      = "#{Odir}/tmp/seq"        # each genome sequence (splitted)
  Ldir      = "#{Odir}/tmp/seq_list"   # contig list
  Idir3k    = "#{Odir}/tmp/seq3k"      # contigs >=3kb
  Idir10k   = "#{Odir}/tmp/seq10k"     # contigs <=10kb
  Rdir      = "#{Odir}/result"         # "#{Odir}/result/refpkg/#{pkg[:name]}/{all,each}/{query,alignment,placement,assign,extract}/#{que[:name]}"
  Fdir      = "#{Odir}/result/fasta"   # "#{Mdir}/C_fixation/CBB/rbcL"
  Logdir    = "#{Odir}/log/tasks"
  HhDb      = $valTools["hhsuite_pfam_db_prefix"]
  CatNrDb   = $valTools["cat_bat_nr_db"]
  CatTxDb   = $valTools["cat_bat_tx_db"]
  VsData    = $valTools["virsorter_data"]   ## Data directory of virsorter

  ## variables
  $fins     = []

  ## Odir exist?
  $stderr.puts "\n\n#{Warmsg} output directory #{Odir} already exists. Overwrite it.\n\n" if OdirExist == "true"

	### run
	NumStep  = Tasks.size
	Tasks.each.with_index(1){ |task, idx|
		Rake::Task[task].invoke(idx)
	}
end
# }}} default (run all tasks)


# {{{ tasks
desc "01-1a.validate_input"
task "01-1a.validate_input", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  ## validate Category
  Cats = %w|taxonomy coverage tetranuc virsorter terL circular| ## valid categories
  cats = Category.gsub(/\s+/, "").gsub("'", "").split(",")
  if cats == %w|all|
    $stderr.puts "\n#{Catemsg} All categories are used to remove contigs. i.e., all of taxonomy,coverage,tetranuc,virsorter,terL,circular\n\n"
  else
    cats.each{ |cat|
      if cat == "all" ## all is used with other categories
        raise("'--category all' should not be combined with other categories. Use only 'all'.")
      end
      unless Cats.include?(cat)
        raise("'--category #{Category}' unknown category '#{cat}'.")
      end
    }
    $stderr.puts "\n#{Catemsg} '#{Category}' are used to remove contigs.\n\n"
  end

  ### validate query and make a copy
  fins = Fins.split(",").inject([]){ |a, path| a += Dir[path.gsub("~", ENV["HOME"])].sort }
  fins.map{ |fin|
    if File.zero?(fin)
      $stderr.puts "#{Warmsg} #{fin} is empty. Skip it."
      next nil
    else
      File.open(fin){ |fr|
        if fr.gets[0] != ">"
          $stderr.puts "#{Warmsg} #{fin} is not valid fasta format. The first character should be '>'."
          next nil
        end
      } 
    end
    fin ## valid fasta
  }.compact

  $stderr.puts ["", "", "\e[1;32m===== check query file (N=#{fins.size}) \e[0m"]
  raise("#{Errmsg} no query file detected.") if fins.size == 0

  idx = 0
  names = {}
  mkdir_p Idir unless File.directory?(Idir)
  $fins  = fins.inject([]){ |a, fin|
    idx += 1
    name = File.basename(fin).split(".")[0..-2]*"."
    raise("#{Errmsg} file name #{name} is not unique.") if names[name]
    names[name] = 1

    ### make copy of input fasta
    path    = "#{Idir}/#{name}.fa"
    path3k  = "#{Idir3k}/#{name}.fa"  ### >=3kb
    path10k = "#{Idir10k}/#{name}.fa" ### <=10kb

    ### make copy of input fasta and extract ge3kb contigs (for virsorter)
    mkdir_p Idir3k unless File.directory?(Idir3k)

    open(path, "w"){ |fw|
      open(path3k, "w"){ |fw3k|
        open(File.absolute_path(fin)){ |fr|
          lab, seq = "", ""
          while l = fr.gets
            if l =~ /^>(\S+)/
              if lab != ""
                seq = seq.gsub(/\s+/, "")
                len = seq.size
                fw.puts ">#{lab}\n#{seq}"
                fw3k.puts ">#{lab}\n#{seq}" if len >= 3000 ### >=3kb
              end

              seq = ""
              lab = l[/^>(\S+)/, 1]
            else
              seq += l.gsub(/\s+/, "")
            end
          end
          ### last entry
          if lab != ""
            seq = seq.gsub(/\s+/, "")
            len = seq.size
            fw.puts ">#{lab}\n#{seq}"
            fw3k.puts ">#{lab}\n#{seq}" if len >= 3000 ### >=3kb
          end
        }
      }
    }

    ### store bp of each fasta (If bp < 20,000, prodigal should use '-p meta'. Otherwise, run railed.)
    bp_path, bp_path3k = 0, 0
    [path, path3k].zip([bp_path, bp_path3k]){ |pa, bp|
      if File.exist?(pa)
        open(pa){ |fr|
          while l = fr.gets
            next if l.strip =~ /^>/
            bp += l.strip.size
          end
        }
      end
    }

    h = { idx: idx, name: name, bp_path: bp_path, bp_path3k: bp_path3k, path: path, path3k: path3k, original: File.absolute_path(fin), checkm_idx: (idx - 1) / ChunkCheckm + 1 }
    $stderr.puts h.inspect

    a << h ### [!] this line must be at the end.
  }
end
desc "01-1b.validate_dbs"
task "01-1b.validate_dbs", ["step"] do |t, args|
  ### hhsuite_db
  %w|a3m.ffdata a3m.ffindex cs219.ffdata cs219.ffindex hhm.ffdata hhm.ffindex|.each{ |suf|
    if Dir["#{HhDb}_#{suf}"].size == 0
      raise("File '#{HhDb}_#{suf}' is not found (hhsuite_pfam_db_prefix: #{HhDb}). Please check config file.")
    end
  }
  if Dir["#{CatNrDb}/*.nr.dmnd"].size == 0
    raise("File '#{CatNrDb}/*.nr.dmnd' is not found (cat_bat_nr_db: #{CatNrDb}). Please check config file.")
  end
  if Dir["#{CatNrDb}/*.nr.fastaid2LCAtaxid"].size == 0
    raise("File '#{CatNrDb}/*.nr.fastaid2LCAtaxid' is not found (cat_bat_nr_db: #{CatNrDb}). Please check config file.")
  end
  if Dir["#{CatTxDb}/names.dmp"].size == 0
    raise("File '#{CatTxDb}/names.dmp' is not found (cat_bat_tx_db: #{CatTxDb}). Please check config file.")
  end
  if Dir["#{CatTxDb}/nodes.dmp"].size == 0
    raise("File '#{CatTxDb}/nodes.dmp' is not found (cat_bat_tx_db: #{CatTxDb}). Please check config file.")
  end

  if Dir["#{CatTxDb}/MAGRE/parsed2.txt"].size == 0
    $stdout.puts "\n\e[1;35m[Creating...]\e[0m '#{CatTxDb}/MAGRE/parsed2.txt'\n"
    sh "ruby #{Sdir}/lib/make_parsed2.rb #{CatTxDb}"
  end
end
desc "01-1c.make_ctg_list"
task "01-1c.make_ctg_list", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  outs = []

	script = "#{Sdir}/#{t.name}.rb"

  mkdir_p Ldir unless File.directory?(Ldir)
  $fins.each{ |fin|
    ### calculate genome size
    fctg  = "#{Ldir}/#{fin[:name]}.ctg.list"

    outs << "ruby #{script} #{fin[:name]} #{fin[:path]} #{fctg}" 
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-2a.tetranucleotide"
task "01-2a.tetranucleotide", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  outs = []

	script = "#{Sdir}/#{t.name}.rb"

  $fins.each{ |fin|
    %w|tetranucleotide|.each{ |_|
      odir = "#{Tdir}/tetranucleotide/#{fin[:name]}"; mkdir_p odir unless File.exist?(odir)
      fout = "#{odir}/tetranucleotide.tsv"
      flog = "#{odir}/tetranucleotide.log"
      next if File.exist?(fout)

      outs << "ruby #{script} #{fin[:path]} #{fout} >#{flog} 2>&1"
    }
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-2b.prepare_for_pca"
task "01-2b.prepare_for_pca", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  outs = []

	script = "#{Sdir}/#{t.name}.rb"

  if File.exist?(Fcov)
    $fins.each{ |fin|
      odir = "#{Tdir}/coverage/#{fin[:name]}"; mkdir_p odir unless File.directory?(odir)
      fout = "#{odir}/coverage.tsv"
      flog = "#{odir}/coverage.log"
      next if File.exist?(fout)
      outs << "ruby #{script} #{fin[:path]} #{Fcov} #{fout} >#{flog} 2>&1"
    }
  else
    $stdout.puts "[!] coverage file is not given. Skip this process."
  end

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-2c.pca"
task "01-2c.pca", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  outs = []

	script = "#{Sdir}/#{t.name}.R"

  $fins.each{ |fin|
    %w|tetranucleotide coverage|.each{ |type|
      odir = "#{Tdir}/#{type}/#{fin[:name]}"
      flog = "#{odir}/#{type}.pca.log"
      fctg = "#{Ldir}/#{fin[:name]}.ctg.list"
      next unless File.exist?("#{odir}/#{type}.tsv") ## fin
      next if     File.exist?("#{odir}/#{type}.pc1_zscore") ## ran?

      ### if just one contig exists, Rscript cause error.
      if IO.readlines(fctg).size == 3 ### header, bin, contig1
        sh "touch #{odir}/#{type}.pc1_zscore" ## create empty file
        next
      end

      outs << "Rscript #{script} #{odir} #{type} >#{flog} 2>&1"
    }
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-3a.prodigal"
task "01-3a.prodigal", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  outs = []

  $fins.each{ |fin|
    ### [a] prodigal
    Codons.each{ |codon| ### condon table: 4, 11
      cleaner = "#{Sdir}/lib/prodigal_cleaner.rb"

      gff  = "cds.gff"
      faa  = "cds.faa"
      fna  = "cds.fna"
      log  = "prodigal.log"
      faa2 = "cds.clean.faa"

      ### [2021-08-27 changed] prodigal10k is built from prodigal
      %w|prodigal|.zip([fin[:path]], [fin[:bp_path]]){ |type, _fin, _bp|
        odir = "#{Tdir}/#{type}/#{fin[:name]}/table#{codon}"; mkdir_p odir unless File.directory?(odir)
        if File.zero?(_fin)
          sh "(cd #{odir} ; touch #{gff} #{faa} #{fna} #{faa2} ; echo 'Input fasta is empty. Empty output files are created.')"; next
        end
        next if File.exist?("#{odir}/#{faa2}")

        if _bp < 20000 ## workaround for "Error:  Sequence must be 20000 characters (only 5824 read)."
          cmd = "prodigal -p meta   -g #{codon} -i #{File.absolute_path(_fin)} -f gff -o #{gff} -a #{faa} -d #{fna} >#{log} 2>&1"
        else
          cmd = "prodigal -p single -g #{codon} -i #{File.absolute_path(_fin)} -f gff -o #{gff} -a #{faa} -d #{fna} >#{log} 2>&1"
        end

        outs << "(cd #{odir} ; #{cmd} && #{cleaner} #{faa} >#{faa2})" 
      }
    }
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-3b.identify_codon_table"
task "01-3b.identify_codon_table", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []
	script = "#{Sdir}/#{t.name}.rb"

  ### SELECT CODON TABLE (either 4 or 11), following checkm
  ## ( ref: https://github.com/Ecogenomics/CheckM/wiki/Genome-Quality-Commands#tree )
  ## Note: The following heuristic is used to establish the translation table used by Prodigal:
  ## use table 11 unless the coding density using table 4 is 5% higher than when using table 11 and the coding density under table 4 is >70%.
  ## Distinguishing between tables 4 and 25 is challenging so CheckM does not attempt to distinguish between these two tables.
  ## If you know the correct translation table for your genomes, it is recommended that you call genes outside of CheckM and
  ## provide CheckM with the protein sequences (see --genes).

  $fins.each{ |fin|
    ### [1] parse prodigal "#{Tdir}/prodigal/#{fin[:name]}/table#{codon}/cds.gff" 
    ### [2] make "#{Tdir}/prodigal/#{fin[:name]}/stats.txt"
    ### [3] make "#{Tdir}/prodigal/#{fin[:name]}/selected" as symlink of either table4 or table11
    ### [4] make "#{Tdir}/prodigal10k/#{fin[:name]}/selected" as symlink of either table4 or table11
    fctg = "#{Tdir}/seq_list/#{fin[:name]}.ctg.list"
    pdir = "#{Tdir}/prodigal/#{fin[:name]}"
    fout = "#{Tdir}/prodigal/#{fin[:name]}/stats.txt"
    flog = "#{Tdir}/prodigal/#{fin[:name]}/stats.log"
    odir1 = "#{Tdir}/prodigal/#{fin[:name]}/selected"
    odir2 = "#{Tdir}/prodigal10k/#{fin[:name]}/selected"
    next if File.directory?(odir1) and File.directory?(odir2) and File.exist?(fout)

    outs << "ruby #{script} #{fctg} #{pdir} #{Codons*","} >#{flog} 2>&1"
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-4.other_tasks"
task "01-4.other_tasks", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []
  signalps = []

  $fins.each{ |fin|
    ### [0] register faa into #$ins
    fin[:faa]    = "#{Tdir}/prodigal/#{fin[:name]}/selected/cds.clean.faa"
    fin[:faa10k] = "#{Tdir}/prodigal10k/#{fin[:name]}/selected/cds.clean.faa"
  }

  ### Constants for pipeline_for_high_sensitive_domain_search
  E_thre  = "1e-10"
  Nreport = 5
  Iter    = 5

  $fins.each{ |fin|
    ### [a] ccfind
    %w|ccfind|.each{ |_|
      script = "#{File.dirname(__FILE__)}/ccfind/ccfind"; raise unless File.exist?(script)

      odir = "#{Tdir}/ccfind/#{fin[:name]}"
      next if File.exist?("#{odir}/result/circ.detected.list") ## ran?
      outs << "#{script} --overwrite --ncpus #{Ncat} #{fin[:path]} #{odir} >/dev/null"
    }

    ### [b] CAT/BAT
    %w|cat|.each{ |_|
      odir = "#{Tdir}/cat_bat/#{fin[:name]}"
      tdir = "#{odir}/tmp"; mkdir_p tdir unless File.directory?(tdir)

      next if File.exist?("#{odir}/BAT.out.bin2classification.txt") and File.exist?("#{odir}/CAT.out.contig2classification.txt") ## ran?

      cmds  = []
      cmds << "CAT bin --force -n #{Ncat} --tmpdir #{tdir} -b #{fin[:path]} -d #{CatNrDb} -t #{CatTxDb} -p #{fin[:faa]} -o #{odir}/BAT.out >/dev/null" ## log file: BAT.out.log
      cmds << "CAT contigs -n #{Ncat} --tmpdir #{tdir} -c #{fin[:path]} -d #{CatNrDb} -t #{CatTxDb} -p #{fin[:faa]} -a #{odir}/BAT.out.alignment.diamond -o #{odir}/CAT.out >/dev/null" ## log file: CAT.out.log
      outs << cmds.join(" && ").gsub(/\s+/, " ")
    }

    ### [c] virsorter
    %w|virsorter|.each{ |_|
      odir = "#{Tdir}/virsorter-ge3kb/#{fin[:name]}"; mkdir_p odir unless File.directory?(odir)
      flog = "#{odir}/virsorter-ge3kb.log"

      next if File.zero?(fin[:path3k])
      next if File.exist?("#{odir}/VIRSorter_global-phage-signal.csv") ## finished?

      outs << "wrapper_phage_contigs_sorter_iPlant.pl -f #{fin[:path3k]} --db 1 --data-dir #{VsData} --wdir #{odir} --ncpu #{Ncat} >#{flog} 2>&1"
    }

    ### [d] pipeline_for_high_sensitive_domain_search
    %w|pipeline_for_high_sensitive_domain_search|.each{ |_|
      script = "#{File.dirname(__FILE__)}/pipeline_for_high_sensitive_domain_search/pipeline_for_high_sensitive_domain_search"; raise unless File.exist?(script)

      odir = "#{Tdir}/pipeline_for_high_sensitive_domain_search/#{fin[:name]}"
      if File.zero?(fin[:faa10k]) ## no <=10kb contig
        mkdir_p  "#{odir}/result/cds.clean.faa/ungrouped" unless File.directory?(odir)
        sh "touch #{odir}/result/cds.clean.faa/ungrouped/besthit.tsv ; echo 'An empty file of #{odir}/result/cds.clean.faa/ungrouped/besthit.tsv is created, because input fasta is empty.'"
        next
      end
      next if File.exist?("#{odir}/result/cds.clean.faa/ungrouped/besthit.tsv") ## ran?

      jdb  = "#{Ddir}/EVG_env_virus.faa" ## proteins included in EVG (Nishimura_et_al_mSphere, 2017), but limited to proteins of >=300 aa
      hhdb = $valTools["hhsuite_pfam_db_prefix"]

      outs << "#{script} --ncpus #{Ncat} --nreport #{Nreport} --iter #{Iter} -i #{fin[:faa10k]} -o #{odir} -J #{jdb} -H #{hhdb} --overwrite >/dev/null"
    }
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, NcatPara, Logdir)
end
desc "01-5.make_contig_table"
task "01-5.make_contig_table", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []
	script = "#{Sdir}/#{t.name}.rb"

  $fins.each{ |fin|
    odir = "#{Rdir}/#{fin[:name]}"; mkdir_p odir unless File.directory?(odir)
    fout = "#{odir}/contig_table.tsv"
    flog = "#{odir}/contig_table.log"
    ftax = "#{CatTxDb}/MAGRE/parsed2.txt"
    outs << "ruby #{script} '#{Category}' #{MinLen} #{Tdir} #{fin[:name]} #{ftax} #{fout} >#{flog} 2>&1"
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
desc "01-6.make_fasta"
task "01-6.make_fasta", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []
	script = "#{Sdir}/#{t.name}.rb"

  $fins.each{ |fin|
    odir = "#{Rdir}/#{fin[:name]}"; mkdir_p odir unless File.directory?(odir)
    fout = "#{odir}/refined.fa"
    flog = "#{odir}/refined.log"
    flst = "#{odir}/contig_table.tsv"
    outs << "ruby #{script} #{fin[:path]} #{flst} #{fout} >#{flog} 2>&1"
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}
