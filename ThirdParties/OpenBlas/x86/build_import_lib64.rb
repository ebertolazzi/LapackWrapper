name = "libopenblas";
cmd1 = "dlltool -z #{name}_x64.def --export-all-symbol x64\\dll\\#{name}_x64.dll"
puts "cmd1 = #{cmd1}\n\n"
system(cmd1)
cmd2 = "lib /machine:X64 /def:#{name}_x64.def"
puts "cmd2 = #{cmd2}\n\n"
system(cmd2)
cmd3 = "move #{name}_x64.lib x64\\lib\\#{name}_x64.lib"
puts "cmd3 = #{cmd3}\n\n"
system(cmd3)
puts "done\n\n"
