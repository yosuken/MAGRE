
fin, fout, *excess = ARGV
raise("argument is not enough") unless fout
raise("argument is too much") if excess.size > 0

I = 4
N = %w|A T G C|

def tetra2idx
  a = [""]
  I.times do 
    b = []
    a.each{ |i| 
      N.each{ |j|
        b << (i + j)
      }
    }
    a = b
  end

  h = {}
  a.each.with_index{ |i, idx| h[i] = idx }
  return h
end

t2i = tetra2idx()
header = %w|sequence| + t2i.keys

open(fout, "w"){ |fw|
  fw.puts header*"\t"

  IO.read(fin).split(/^>/)[1..-1].each{ |ent|
    cnt = Array.new(t2i.size, 0)

    lab, *seq = ent.split("\n")
    lab = lab.split(/\s+/)[0]
    seq = seq.join("").gsub(/\s+/, "").upcase
    m   = seq.size - I + 1 ## begin postion of the last frame

    (0..m).each{ |i|
      f = seq[i, I]
      idx = t2i[f] 
      cnt[idx] += 1 if idx
    }

    fw.puts [lab, cnt]*"\t"
  }
}
