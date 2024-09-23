name = "openblas";
cmd1 = "dlltool -z #{name}_x86.def --export-all-symbol dll\\#{name}_x86.dll"
puts "cmd1 = #{cmd1}\n\n"
system(cmd1)
cmd2 = "lib /machine:x86 /def:#{name}_x86.def"
puts "cmd2 = #{cmd2}\n\n"
system(cmd2)
cmd3 = "copy #{name}_x86.lib lib\\#{name}_x86.lib"
puts "cmd3 = #{cmd3}\n\n"
system(cmd3)
puts "done\n\n"
