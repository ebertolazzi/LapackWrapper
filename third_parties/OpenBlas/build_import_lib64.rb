name = ARGV[0];
cmd1 = "dlltool -z #{name}_x64.def --export-all-symbol #{name}.dll"
cmd2 = "lib /machine:X64 /def:#{name}_x64.def"
puts "cmd1 = #{cmd1}\n\n"
system(cmd1)
puts "cmd2 = #{cmd2}\n\n"
system(cmd2)
puts "done\n\n"
