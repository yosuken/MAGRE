#!/usr/bin/env ruby

fin, *excess = ARGV
raise("argument is not enough") unless fin
raise("argument is too much") if excess.size > 0

## [1] trim sequential "X" or "*" at the C-terminal and
## [2] remove "\n" in sequence

IO.read(fin).split(/^>/)[1..-1].each{ |ent|
  lab, *seq = ent.split("\n")
  seq = seq.join("").gsub(/[X\*]+$/, "").gsub(/\s+/, "")
  puts [">#{lab}", seq]
}
