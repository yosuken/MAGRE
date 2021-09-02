
category, minlen, Tdir, bin, ftax, fout, *excess = ARGV 
raise("argument is not enough") unless fout
raise("argument is too much") if excess.size > 0

### parse category
cats = %w|taxonomy coverage tetranuc virsorter terL circular|
category = category.gsub(/\s+/, "").gsub("'", "").split(",")
category = category == ["all"] ? cats : category
Category = category + ["length"] ### add length filter
MinLen   = minlen.to_i


# {{{ setting and documentation
###                   ###
### [1] Cutoff values ###
###                   ###
MinCovPC1  = -2.5  ## PC1 of coverage
MaxCovPC1  =  2.5  ## PC1 of coverage
MinTetPC1  = -2.5  ## PC1 of tetranucleotide
MaxTetPC1  =  2.5  ## PC1 of tetranucleotide
MaxVirsFrc =  0.5  ## max virsorter detected viral fraction to recognize as cellular
MinTerLlen = 10000 ## min length of contig to recognize as cellular by terL detection
TerL_prob  = 97    ## probability threshold of hhsearch
### Pfam terL model (including flagment model)
TerLs = %w|Terminase_1 Terminase_3 Terminase_6 Terminase_GpA DNA_pack_N Terminase_3C Terminase_6C|.inject({}){|h,i|h[i]=1;h}

### make params log
puts "MinCovPC1: #{MinCovPC1}"
puts "MaxCovPC1: #{MaxCovPC1}"
puts "MinTetPC1: #{MinTetPC1}"
puts "MaxTetPC1: #{MaxTetPC1}"
puts "MaxVirsFrc: #{MaxVirsFrc}"
puts "MinTerLlen: #{MinTerLlen}"
puts "TerL_prob: #{TerL_prob}"
puts "TerLs: #{TerLs.keys*", "}"
$stdout.flush

###                                     ###
### [2] Taxonomy refinement (polishing) ###
###                                     ###
# (0) trim starred taxa 
# Starred taxa should be treated carefully. (ref: https://github.com/dutilh/CAT#marking-suggestive-classifications-with-an-asterisk)
#
# (1) trim species and lower annotation 
# Spiecies annotation is too stringent, becauses species A (from MAG) and species B (from isolate or SAG or MAG) might belong to same species.
# In addition, unspecific species annotation (e.g., Alphaproteobacteria bacterium [species] and Candidatus Poseidoniales archaeon [species]) should be removed.
#
# (2) trim "unclassified xxx" (e.g., unclassified Pelagibacteraceae) and "environmental samples"  recursively from a leaf to the root. 
# When non-unclassified taxonomy is leaf, stop this trimming.

###                          ###
### [3] taxonomic filtering  ###
###                          ###
# If (refined) taxon of contig is NOT
#   - ancestor of bin taxon
#   and
#   - descendant of bin taxon
# then remove.
#
# [example]
# bin:  A;B;C;D
# ctg1: A;B;C      --> keep
# ctg2: A;B;C;D;E  --> keep
# ctg3: A;B;C;F    --> remove
# ctg4: A;G;H      --> remove
# }}}


# {{{ files
fctg  = "#{Tdir}/seq_list/#{bin}.ctg.list"
fcds  = "#{Tdir}/prodigal/#{bin}/stats.txt"              ### prodigal (for codon table)
fprod = "#{Tdir}/prodigal/#{bin}/selected/cds.gff"       ### prodigal (coding len, density)
fcirc = "#{Tdir}/ccfind/#{bin}/result/circ.detected.list" ### ccfind circular detected contigs
fvirs = "#{Tdir}/virsorter-ge3kb/#{bin}/Predicted_viral_sequences/VIRSorter_*-?.gb" ### virsorter genbank file
fvirc = "#{Tdir}/virsorter-ge3kb/#{bin}/fasta/VIRSorter_circu.list"                 ### virsorter circular list
fterL = "#{Tdir}/pipeline_for_high_sensitive_domain_search/#{bin}/result/cds.clean.faa/ungrouped/besthit.tsv" ### terL detected
fcat  = "#{Tdir}/cat_bat/#{bin}/CAT.out.contig2classification.txt"
fbat  = "#{Tdir}/cat_bat/#{bin}/BAT.out.bin2classification.txt"
fcov  = "#{Tdir}/coverage/#{bin}/coverage.pc1_zscore"
ftet  = "#{Tdir}/tetranucleotide/#{bin}/tetranucleotide.pc1_zscore"
# }}}


# {{{ functions
class Range
  include Comparable

  def <=>(other)
    self.min <=> other.min
  end
  def overlap?(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return false
    else
      return true
    end
  end
  def &(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return nil
    else
      return ([s.min, t.min].max)..([s.max, t.max].min)
    end
  end
  def |(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return [s, t]
    else
      return [([s.min, t.min].min)..([s.max, t.max].max)]
    end
  end
  def include?(t)
    s = self
    o = s & t  ## overlap of s and t
    return false unless o
    if o.min == t.min and o.max == t.max ## s includes t
      return true
    else
      return false
    end
  end
end
def merge_ranges(ranges) ## ranges = [3..300, 320..500, 504..732, ...]
  return ranges if ranges.size < 2

  rs  = ranges.sort
  (-rs.size..-2).each{ |j| ## index from right
    merged = rs[j] | rs[j+1]
    case merged.size
    when 1 ## overlap detected
      rs[j] = merged[0]
      rs.delete_at(j+1)
    when 2 ## overlap not detected
      ## do nothing
    else raise
    end
  }
  return rs
end
def fill_blank(ctg2info, posi)
  ctg2info.each{ |ctg, info|
    (posi - info.size).times do info << "-" end
    (0..info.size-1).each{ |i| info[i] ||= "-" }
  }
end
def tax_refine(a)
  # a: %w|A B C D|
  ## (1) trim species and lower taxon
  flag = 0
  a = a.map{ |i|
    flag = 1 if i =~/ \[species\]$/
    flag == 1 ? nil : i ## discard after species
  }
  a = a.compact
  a = a.reverse

  ## (2) trim "unclassified xxx" recursively from a leaf to the root.
  flag = 0
  a = a.map{ |i|
    flag = 1 if i !~ /^unclassified / and i !~ /^environmental samples$/ ## break if non-unclassified taxonomy comes.
    flag == 0 ? nil : i
  }
  a = a.compact
  a = a.reverse
  return a
end
# }}} functions


header = %w|bin_or_contig
length  gc  gc_skew  percent_N
translTable  codingLen  codingDens
cat_bat_classification  cat_bat_reason  cat_bat_lineage  cat_bat_lineage_scores  cat_bat_lineage_names  polished_lineage  polished_lineage_names 
coverage_pc1  tetranucleotide_pc1
ccfind_circular
virsorter_category  virsorter_viral_length  virsorter_viral_position  virsorter_contig_length  virsorter_viral_fraction
terL_gene  terL_aa_len  terL_besthit  terL_best_prob
all_categories_pass_threshold  selected_categories  is_contamination
|
ctg2info = {}
virS2ctg = {} ### virsorter contig (PACC02000035_1 --> PACC02000035.1)
ctg2contam = Hash.new{ |h, i| h[i] = [] }
posi = 0


# {{{ list of contig (n=4)
size = 4
%w|list|.each{ |_|
  IO.readlines(fctg)[1..-1].each{ |l|
    ### include bin and all contigs
    # #name length GC GC-Skew A T G C N
    ctg, len, gc, skew, ns = l.chomp.split("\t").values_at(0..3, 8)
    ns = "%.6g" % (100.0 * ns.to_i / len.to_i)
    ctg2info[ctg] = [len, gc, skew, ns]
    virS = ctg.gsub(".", "_")
    virS2ctg[virS] = ctg

    ctg2contam[ctg] << "length:#{len}" if len.to_i < MinLen
  }
}
posi += size
fill_blank(ctg2info, posi)
# }}}


# {{{ prodigal (n=3)
size  = 3
%w|prodigal|.each{ |_|
  table  = IO.read(fcds)[/^selected codon table: (\d+)/, 1]

  ranges = []
  ctg = ""
  # bin_table, bin_codingLen, bin_codingDens = nil, 0, 0
  Dir[fprod].each{ |fin|
    next if File.zero?(fin)
    IO.readlines(fin).each{ |l|
      next if l =~ /^#/
      _ctg, start, stop = l.chomp.split("\t").values_at(0, 3, 4)

      if ctg != _ctg  ### if contig changed
        if ctg != "" ### skip first line
          ranges       = merge_ranges(ranges)
          codingLen    = ranges.map{ |r| r.max - r.min + 1 }.inject(&:+)
          codingDens   = "%.7g" % (codingLen.to_f / ctg2info[ctg][0].to_i)
          ctg2info[ctg][posi]   = table
          ctg2info[ctg][posi+1] = codingLen
          ctg2info[ctg][posi+2] = codingDens
          ### bin stat
          # bin_codingLen += codingLen
        end
        ranges       = []
        ctg           = _ctg
      end

      ranges << (start.to_i..stop.to_i)
    }

    ## process last contig
    if ranges.size > 0 and ctg != ""
      # p [ctg, ranges]
      ranges       = merge_ranges(ranges)
      codingLen    = ranges.map{ |r| r.max - r.min + 1 }.inject(&:+)
      codingDens   = "%.7g" % (codingLen.to_f / ctg2info[ctg][0].to_i)
      ctg2info[ctg][posi]   = table
      ctg2info[ctg][posi+1] = codingLen
      ctg2info[ctg][posi+2] = codingDens

      ### bin stat
      ctg2info[bin][posi]   = table ## bin_table
      ctg2info[bin][posi+1] = ctg2info.map{ |ctg, info| info[posi+1] }.compact.inject(&:+) ### bin_codingLen
      ctg2info[bin][posi+2] = "%.6g" % (ctg2info[bin][posi+1].to_f / ctg2info[bin][0].to_i)  ### bin_codingdens
    end
  }
}
posi += size
fill_blank(ctg2info, posi)
# }}}


# {{{ cat_bat (n=7)
size = 7  ## classification  reason  lineage  lineage_scores  lineage_names  polished_lineage  polished_lineage_names
%w|cat_bat|.each{ |_t|
  ### parse taxonomy (parsed2.txt) [a]
  tid2info = {}
  IO.readlines(ftax).each{ |l|
    ## 2   Bacteria|superkingdom|2   cellular organisms|no rank|131567   root|no rank|1
    tid, info = l.chomp.split("\t")
    name, rank = info.split("|")[0..1] ## remove taxid
    if rank == "no rank" or rank == "superkingdom"
      tid2info[tid] = "#{name}"
    else
      tid2info[tid] = "#{name} [#{rank}]"
    end
  }

  ### parse BAT [b-1]
#              # bin   classification                   reason                lineage            lineage scores
# ERS477931_bin.1.fa       classified  based on 1967/2020 ORFs  1;131567;2;2323;49928  1.00;1.00;1.00;0.97;0.97
  Dir[fbat].each{ |fin|
    next if File.zero?(fin)
    IO.readlines(fin)[1..-1].each{ |l|
      a = l.chomp.split("\t")
      _bin = a[0].split(".")[0..-2]*"." ## remove ".fa"
      raise if _bin != bin

      lin_ids     = a[3]||""
      lin_scrs    = a[4]||""
      lin_ids     = lin_ids.split(";").select{ |i| i =~ /^\d+$/ }*";"        ## remove starred taxa
      lin_names   = lin_ids.split(";").map{ |tid| tid2info[tid] }*" ; "      ## lineage_names
      r_lin_names = tax_refine(lin_names.split(" ; "))*" ; "                 ## taxon refinement
      r_lin_ids   = lin_ids.split(";")[0, r_lin_names.split(" ; ").size]*";" ## same size as r_lin_names

      ctg2info[_bin] += [a[1], a[2], lin_ids, lin_scrs, lin_names, r_lin_ids, r_lin_names]  ## bin entry
    }
  }

  ### parse CAT [b-2]
  #        0               1       2        3               4
  # # contig  classification  reason  lineage  lineage scores
  Dir[fcat].each{ |fin|
    next if File.zero?(fin)
    IO.readlines(fin)[1..-1].each{ |l|
      a = l.chomp.split("\t")
      ctg = a[0]

      lin_ids     = a[3]||""
      lin_scrs    = a[4]||""
      lin_ids     = lin_ids.split(";").select{ |i| i =~ /^\d+$/ }*";"        ## remove starred taxa
      lin_names   = lin_ids.split(";").map{ |tid| tid2info[tid] }*" ; "      ## lineage_names
      r_lin_names = tax_refine(lin_names.split(" ; "))*" ; "                 ## taxon refinement
      r_lin_ids   = lin_ids.split(";")[0, r_lin_names.split(" ; ").size]*";" ## same size as r_lin_names

      ### detect contamination (taxonomic inconsistency)
      b = ctg2info[bin][posi+6] ## polished taxonomy of bin
      c = r_lin_names           ## polished taxonomy of contig
      contam  = (b[0, c.size] == c or c[0, b.size] == b) ? [] : ["taxonomy"]
      ctg2contam[ctg] += contam

      ctg2info[ctg] += [a[1], a[2], lin_ids, lin_scrs, lin_names, r_lin_ids, r_lin_names]  ## contigs
    }
  }
}
posi += size
fill_blank(ctg2info, posi)
# }}}


# {{{ coverage and tetranucleotide (n=2)
size = 2
[fcov, ftet].each.with_index{ |fin, idx|
  next unless File.exist?(fin)
  next if File.zero?(fin)
  IO.readlines(fin)[1..-1].each{ |l|
    ctg, val = l.chomp.split("\t")
    val = val.to_f
    raise("unknown contig: #{ctg} in #{bin}") if ctg2info[ctg].size == 0
    ctg2info[ctg][posi + idx] = "%.7g" % val

    ### contam?
    case idx
    when 0 ### coverage
      ctg2contam[ctg] << "coverage:#{val}" if val < MinCovPC1 or val > MaxCovPC1
    when 1 ### tetranucleotide
      ctg2contam[ctg] << "tetranuc:#{val}" if val < MinTetPC1 or val > MaxTetPC1
    else raise
    end
  }
}
posi += size
fill_blank(ctg2info, posi)
# }}}


# {{{ ccfind (n=1)
size  = 1
%w|ccfind|.each{ |_t|
  Dir[fcirc].each{ |fin|
    next if File.zero?(fin)
    IO.readlines(fcirc).each{ |l|
      ctg = l.chomp.split("\t")[0]
      ctg2info[ctg][posi] = "circular"
      ctg2contam[ctg] << "circular"
    }
  }
}
posi += size
fill_blank(ctg2info, posi)
# }}}


# {{{ virsorter (n=5)
size  = 5
%w|virsorter|.each{ |_t|
  Dir[fvirs].each{ |fin|
    next if File.zero?(fin)
    category = File.basename(fin)[/-(\d)\.gb$/, 1] ##VIRSorter_*-?.gb
    IO.readlines(fin).each{ |l|
      # LOCUS       VIRSorter_SRS992687_N0001491         3028 bp    dna     linear   ENV 12/11/19
      # LOCUS       VIRSorter_SRS992687_N0000009-circular        35568 bp    dna     circular ENV 12/11/19
      # LOCUS       VIRSorter_SRS992687_N0000143_gene_1_gene_18-0-8109-cat_6         8109 bp    dna     linear   ENV 12/11/19
      if l =~ /^LOCUS/
        lab, vlen = l.chomp.split(/\s+/).values_at(1, 2)
        ctg = case lab
             when /^VIRSorter_(\S+)-circular/ then $1 ## circular --> [!!!] length is calculated after removal of terminal redundancy! (could be category 1-6)
             when /^VIRSorter_(\S+)_gene_\d+_gene_\d+-\d+-\d+-cat_\d/ then $1 ## provirus (could be category 4-6)
             when /^VIRSorter_(\S+)$/ then $1 ## virus (could be category 1-3)
             else raise
             end
        po = case lab
             when /^VIRSorter_(\S+)_gene_\d+_gene_\d+-(\d+)-(\d+)-cat_\d/ then "#{$2.to_i+1}-#{$3}" ## provirus (could be category 4-6)
             else "1-#{vlen}" ## virus (could be category 1-3)
             end

        ### [!!!] virsorter contigname workaround
        ### e.g.) PACC02000035_1 --> PACC02000035.1
        ### currently, only the modification ("." -> "_") is rescued. Other conversion is missing if exists.
        unless ctg2info[ctg]
          ctg = virS2ctg[ctg]
          unless ctg2info[ctg]
            p [ctg, l]
            raise("#{ctg}: not valid contig name.")
          end
        end

        unless ctg2info[ctg][posi] 
          ctg2info[ctg][posi]   = category
          ctg2info[ctg][posi+1] = vlen
          ctg2info[ctg][posi+2] = po
          ctg2info[ctg][posi+3] = ctg2info[ctg][0] ## contig length
        else ## not 1st part --> must be provirus
          ctg2info[ctg][posi]   += ", #{category}" ## add
          ctg2info[ctg][posi+1] += ", #{vlen}"     ## add
          ctg2info[ctg][posi+2] += ", #{po}"       ## add
        end
      end
    }
  }
  Dir[fvirc].each{ |fin|
    next if File.zero?(fin)
    IO.readlines(fin).each{ |l|
      ## VIRSorter_ERS477931_N0000002-circular   188442
      lab, len = l.split(/\s+/)
      ctg = lab[/^VIRSorter_(\S+)-circular$/, 1]
      unless ctg2info[ctg]
        ctg = virS2ctg[ctg]
        unless ctg2info[ctg]
          p [ctg, l]
          raise("#{ctg}: not valid contig name.")
        end
      end
      ctg2info[ctg][posi+3] = len ## contig length (without terminal redundancy)
    }
  }

  ### fraction of viral part
  ctg2info.each{ |ctg, info|
    if info[posi] and info[posi] != "-"
      val = info[posi+1].split(", ").map(&:to_f).inject(&:+) / info[posi+3].to_f ## fraction of viral part
      ctg2info[ctg][posi+4] = "%.4f" % val
      ctg2contam[ctg] << "virsorter:#{"%.4f" % val}" if val > MaxVirsFrc
    end
  }
}
posi += size
fill_blank(ctg2info, posi)
# }}}


# {{{ terL (n=4)
size  = 4
%w|terL|.each{ |_t|
  Dir[fterL].each{ |fin|
    next if File.zero?(fin)
    IO.readlines(fin)[1..-1].each{ |l|
      pid, plen, hmm_name, prob = l.chomp.split("\t").values_at(0, 1, 6, 8)
      next unless TerLs[hmm_name]    ### terL motif (including flagment motif)
      next if prob.to_f <= TerL_prob ### probability larger than threshold
      ctg = pid.split("_")[0..-2]*"_"

      unless ctg2info[ctg][posi] 
        ctg2info[ctg][posi]   = pid
        ctg2info[ctg][posi+1] = plen
        ctg2info[ctg][posi+2] = hmm_name
        ctg2info[ctg][posi+3] = prob
      else
        ctg2info[ctg][posi]   += ", #{pid}"      ## add
        ctg2info[ctg][posi+1] += ", #{plen}"     ## add
        ctg2info[ctg][posi+2] += ", #{hmm_name}" ## add
        ctg2info[ctg][posi+3] += ", #{prob}"     ## add
      end
    }
  }

  ## if contig length is less than 10000, terL detected contig is regarded as viral.
  ctg2info.each{ |ctg, info|
    if info[0].to_i < MinTerLlen and info[posi] and info[posi] != "-"
      ctg2contam[ctg] << "terL"
    end
  }
}
posi += size
fill_blank(ctg2info, posi)
# }}}


# {{{ contamination (n=3)
size  = 3
p Category
ctg2contam.each{ |ctg, contam|
  if contam.size > 0
    _contam = [] ## selected contam
    contam.each{ |c|
      _contam << c if Category.include?(c.split(":")[0]) ### if selected category
    }
    if _contam.size > 0 ## if selected category is not zero
      ctg2info[ctg] += [contam*";", _contam*";", "Y"]
    else
      ctg2info[ctg] += [contam*";", "-", "-"]
    end
  else
    ctg2info[ctg] += %w|- - -|
  end
}
posi += size
fill_blank(ctg2info, posi)
# }}}


### make output
open(fout, "w"){ |fw|
  fw.puts header*"\t"
  ctg2info.each{ |ctg, a| fw.puts [ctg, a]*"\t" }
}
