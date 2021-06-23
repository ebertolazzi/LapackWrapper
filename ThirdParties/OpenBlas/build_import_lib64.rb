name = "libopenblas";
#cmd1 = "dlltool -z #{name}_x64.def --export-all-symbol #{name}.dll"
#puts "cmd1 = #{cmd1}\n\n"
#system(cmd1)
cmd2 = "lib /machine:X64 /def:#{name}_x64.def"
puts "cmd2 = #{cmd2}\n\n"
system(cmd2)
puts "done\n\n"
