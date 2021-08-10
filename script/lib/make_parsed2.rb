#!/usr/bin/env ruby

require 'rake'

usage = <<EOS
[usage]
ruby script/02.make_parsed2.rb <idir>

[example]
ruby script/02.make_parsed2.rb Data/2016/20160512/taxdump

[requirement]
1) <idir> should be directory.
2) <idir>/names.dmp should exist.
2) <idir>/nodes.dmp should exist.
EOS

idir = ARGV[0]
(puts usage; exit) if ARGV.size != 1

fin1  = "#{idir}/names.dmp"
fin2  = "#{idir}/nodes.dmp"
(puts usage; exit) unless File.exist?(fin1)
(puts usage; exit) unless File.exist?(fin2)

odir  = "#{idir}/mag_refine"; mkdir_p odir unless File.directory?(odir)
fout  = open("#{odir}/parsed2.txt", "w")
ranks = {} # ranks["family"]["31"] = 1

id2parent = {}
id2name   = {}
id2rank   = {}

IO.readlines(fin1).each{ |l|
	id, name, type = l.chomp.split(/\t/).values_at(0, 2, 6)
	next if type != "scientific name"
	#name_ary += [id, name]
	id2name[id] = name
}
puts "parsed: name file"

IO.readlines(fin2).each{ |l|
	id, parent_id, rank = l.chomp.split(/\t/).values_at(0, 2, 4)
	raise if id2parent[id] # overlap check
	id2parent[id] = parent_id #if id != "1" # "1" => "1"
	id2rank[id]   = rank
}
puts "parsed: node file"

id2parent.each_key{ |input|
	output = {} # { "kingdom_id" => "2", "kingdom_name" => "Bacteria", ...}
	id = input.dup
	idx = 0
	ids = []
	while parent = id2parent[id]
		ids << id
    break if id == "1" and parent == "1"
		id = parent
	end
	fout.puts [input, ids.map{ |id| [id2name[id], id2rank[id], id]*"|" }*"\t"]*"\t"
}
