
fin, flst, fout, *excess = ARGV 
raise("argument is not enough") unless fout
raise("argument is too much") if excess.size > 0

remove = {}
IO.readlines(flst)[2..-1].each{ |l|
  ctg, flag = l.chomp.split("\t").values_at(0, 29)
  remove[ctg] = 1 if flag == "Y"
}

s1, c1, s2, c2 = 0, 0, 0, 0 
open(fout, "w"){ |fw|
  IO.read(fin).split(/^>/)[1..-1].each{ |ent|
    lab, *seq = ent.split("\n")
    ctg = lab[/^(\S+)/, 1]
    len = seq.join.gsub(/\s+/, "").size

    s1 += len
    c1 += 1

    next if remove[ctg] ## skip if contamination flagged contig

    fw.puts ">#{ent}"
    s2 += len
    c2 += 1
  }
}

### write statistics to log file
sp = "%.4g" % (s2 * 100.0 / s1)
cp = "%.4g" % (c2 * 100.0 / c1)
puts "genome size: #{s1} --> #{s2} (#{sp}% retained)"
puts "number of contigs: #{c1} --> #{c2} (#{cp}% retained)"
