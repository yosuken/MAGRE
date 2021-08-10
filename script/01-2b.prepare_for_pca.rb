
fin, Fcov, fout, *excess = ARGV
raise("argument is not enough") unless fout
raise("argument is too much") if excess.size > 0

ctgs = {}
IO.read(fin).split(/^>/)[1..-1].each{ |ent|
  lab, *seq = ent.split("\n")
  ctg = lab.split(/\s+/)[0]
  ctgs[ctg] = 1
}

type = "" ### jgi or tsv
out  = []
size = -1
IO.readlines(Fcov).each.with_index{ |l, idx|
  if idx == 0
    if l =~ /^contigName\tcontigLen\ttotalAvgDepth\t\S+.bam\t\S+.bam-var/
      ### format of "jgi_summarize_bam_contig_depths"
      ### contigName      contigLen       totalAvgDepth   ERR315856.bam   ERR315856.bam-var
      ### --> cut -f1,4,6,8,...
      type = "jgi"
    else ### regard as TSV of "contig", "cov1", "cov2", ...
      type = "tsv"
    end
  end

  a = l.chomp.split("\t")
  if type == "jgi"
    n_smpl = (a.size - 3) / 2
    idxs = (0...n_smpl).map{ |i| i * 2 + 3 }
    a = a.values_at(0, *idxs)
  end

  if idx == 0
    ### [!!!] header is REQUIRED
    a[0] = "sequence" ### 1st column of header
    size = a.size     ### register tsv size
  else
    next unless ctgs[a[0]] ### ignore if the first column is not a name of contigs of the bin
    ctgs.delete(a[0]) ### remove stored contig name
  end

  if size != a.size ### validate number of columns in TSV
    raise("invalid format: #{Fcov}; line #{idx+1} has #{a.size} columns, but line 1 is #{sizes} columns.")
  end

  out << a*"\t"
}

if ctgs.size > 0 ### not all contigs were detected
  raise("\e[1;31m[!!!] Error:\e[0m NOT ALL contigs were detected in #{Fcov}, for example: #{ctg.first[0]}")
else
  open(fout, "w"){ |fw| fw.puts out }
end


