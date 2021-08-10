
bin, fin, fctg, *excess = ARGV
raise("argument is not enough") unless fctg
raise("argument is too much") if excess.size > 0

ctg2info = {} ## ctg, len, gc, skew, a, t, g, c, n
IO.read(fin).split(/^>/)[1..-1].each{ |ent|
  lab, *seq = ent.split("\n")
  seq = seq.join("").gsub(/\s+/, "").upcase

  ctg = lab.split(/\s+/)[0]
  len = seq.size

  ### nuc counter
  info = [len]
  h = Hash.new(0)
  (0..seq.size-1).each{ |i| h[seq[i]] += 1 }
  info << "%.6g" % ((h["G"] + h["C"]) * 100.0 / len) ## GC
  info << "%.6g" % ((h["G"] - h["C"]) / (h["G"] + h["C"]).to_f) ## GC-skew
  info <<  h["A"] ## A
  info <<  h["T"] ## T
  info <<  h["G"] ## G
  info <<  h["C"] ## C
  info <<  h["N"] ## N

  ctg2info[ctg] = info
}

open(fctg, "w"){ |fw|
  fw.puts %w|contig length gc gc_skew A T G C N|*"\t"
  ### all contigs sum
  info = [0, "", "", 0, 0, 0, 0, 0]
  ctg2info.values.each{ |a|
    (0..info.size-1).each{ |i|
      info[i] += a[i]
    }
  }
  info[1] = "%.6g" % ((info[-3] + info[-2]) * 100.0 / info[0]) ## GC
  info[2] = "%.6g" % ((info[-3] - info[-2]) / (info[-3] + info[-2]).to_f) ## GC-skew
  fw.puts [bin, info]*"\t"

  ### each contig
  ctg2info.sort_by{ |ctg, info| [-info[0].to_i, ctg] }.each{ |ctg, info| ## ordered by length (descendant)
    fw.puts [ctg, info]*"\t"
  }
}
