
require 'rake'
Fctg, Pdir, codons, *excess = ARGV
raise("argument is not enough") unless codons
raise("argument is too much") if excess.size > 0

bin  = File.basename(Pdir)
codons = codons.split(",") ## [!!!] codons should be "4,11"
fout = "#{Pdir}/stats.txt"
fgffs = []; codons.each{ |codon| fgffs << "#{Pdir}/table#{codon}/cds.gff" }

def parse_coding_len(fgff)
  c_len = 0
  open(fgff){ |fr|
    while l = fr.gets
      next if l[0] == "#" ### comment line
      s, e = l.chomp.split("\t").values_at(3, 4).map(&:to_i)
      c_len += (e - s + 1)
    end
  }
  return c_len
end
def make_prodigal10k()
  figff = "#{Pdir}/selected/cds.gff"
  fifaa = "#{Pdir}/selected/cds.clean.faa"

  _pdir = Pdir.gsub(%r{/prodigal/}, "/prodigal10k/") ## prodigal result for le10kb contigs
  odir  = "#{_pdir}/selected"; mkdir_p odir unless File.directory?(odir)
  fogff = "#{odir}/cds.gff"
  fofaa = "#{odir}/cds.clean.faa"

  ### parse contig info
  ctg2len = {}
  IO.readlines(Fctg)[2..-1].each{ |l|
    # contig  length  gc      gc_skew A       T       G       C       N
    ctg, len = l.chomp.split("\t").values_at(0, 1)
    ctg2len[ctg] = len.to_i
  }

  ### parse gff
  flag = 0
  open(fogff, "w"){ |fw|
    IO.readlines(figff).each{ |l|
      if l =~ /^##gff/ ## header
        fw.puts l
      elsif l =~ /^# / and l =~ /seqhdr="([^"\s]+)/ ## parse name
        ctg = $1
        if ctg2len[ctg] <= 10000 ### <=10kb
          flag = 1
          fw.puts l
        else
          flag = 0
        end
      else 
        fw.puts l if flag == 1
      end
    }
  }

  ### parse faa
  open(fofaa, "w"){ |fw|
    IO.read(fifaa).split(/^>/)[1..-1].each{ |ent|
      ctg = ent.split("\n")[0].split(/\s+/)[0].split("_")[0..-2]*"_" ## <ctg>_1, <ctg>_2 --> <ctg>
      if ctg2len[ctg] <= 10000 ### <=10kb
        fw.puts ">#{ent}"
      end
    }
  }
end

info = [] ## tot_len, table4_coding, table11_coding, %table4_coding, %table11_coding
out  = []

### [0] total genome size
size = IO.readlines(Fctg)[1].split("\t")[1] ## bin length
info << size.to_i
out << "genome length: #{size}"

### [1] parse prodigal "#{Tdir}/prodigal/#{fin[:name]}/table#{codon}/cds.gff"  (table4, table11)
fgffs.each{ |fgff|
  codon = fgff.split("/")[-2]
  c_len = parse_coding_len(fgff) 
  info << c_len
  out << "#{codon} coding length: #{c_len}"
}

### [2] make           "#{Tdir}/prodigal/#{fin[:name]}/stats.txt"
### [3] make           "#{Tdir}/prodigal/#{fin[:name]}/selected" as symlink of either table4 or table11
info[3] = 100.0 * info[1] / info[0] ## %table4_coding
info[4] = 100.0 * info[2] / info[0] ## %table11_coding
out << "%table4 coding density: #{info[3]}"
out << "%table11 coding density: #{info[4]}"

## ( ref: https://github.com/Ecogenomics/CheckM/wiki/Genome-Quality-Commands#tree )
## Note: The following heuristic is used to establish the translation table used by Prodigal:
## use table 11 unless the coding density using table 4 is 5% higher than when using table 11 and the coding density under table 4 is >70%.
if info[3] > 70 and info[3] - info[4] >= 5
  info[5] = "4"
  out << "selected codon table: 4"
  sh "(cd #{Pdir} ;  ln -s table4 selected)" unless File.directory?("#{Pdir}/selected")

  ### make prodigal10k
  # sh "(cd #{_pdir} ; ln -s table4 selected)" 
  make_prodigal10k("4")
else
  info[5] = "11"
  out << "selected codon table: 11"
  sh "(cd #{Pdir} ;  ln -s table11 selected)" unless File.directory?("#{Pdir}/selected")

  # sh "(cd #{_pdir} ; ln -s table11 selected)" 
end

### make prodigal10k
make_prodigal10k()

### make stats.txt
open(fout, "w"){ |fw| fw.puts out }

