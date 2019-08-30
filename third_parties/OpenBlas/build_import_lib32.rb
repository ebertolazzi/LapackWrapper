name = ARGV[1];
sh "dlltool -z #{name}_x86.def --export-all-symbol #{name}.dll"
sh "lib /machine:X86 /def:#{name}_x86.def"
