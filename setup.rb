puts "Setup submodules"
system("git submodule update --init --recursive --checkout")
system('git submodule foreach --recursive git submodule init')
system('git submodule foreach --recursive git submodule update')
system('git submodule foreach --recursive git submodule sync')
system('(cd cmake_utils; git checkout main)')
system('(cd submodules/UtilsLite; git checkout master)')
system('(cd submodules/UtilsLite/cmake_utils; git checkout main)')

puts ARGV

if ARGV.size() > 0 && ARGV[0] == "--last" then
  puts "\nUpdate submodules to last version"
  system('git submodule foreach --recursive git pull')
end
