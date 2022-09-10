%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

require 'rake/clean'

CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }

case RUBY_PLATFORM
when /darwin/
  OS = :mac
when /linux/
  OS = :linux
when /cygwin|mswin|mingw|bccwin|wince|emx/
  OS = :win
when /msys/
  OS = :mingw
end

require_relative "./Rakefile_common.rb"

file_base = File.expand_path(File.dirname(__FILE__)).to_s

cmd_cmake_build = ""
if COMPILE_EXECUTABLE then
  cmd_cmake_build += ' -DEB_ENABLE_TESTS:VAR=ON '
else
  cmd_cmake_build += ' -DEB_ENABLE_TESTS:VAR=OFF '
end
if COMPILE_DYNAMIC then
  cmd_cmake_build += ' -DEB_BUILD_SHARED:VAR=ON '
else
  cmd_cmake_build += ' -DEB_BUILD_SHARED:VAR=OFF '
end
if COMPILE_DEBUG then
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=STATUS '
else
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=STATUS '
end

desc "default task --> build"
task :default => :build

FileUtils.cp './cmake/CMakeLists-cflags.txt', 'submodules/Utils/cmake/CMakeLists-cflags.txt'

task :default => [:build]

desc "build LAPACK_WRAPPER"
task :build do
  case OS
  when :mac
    puts "LAPACK_WRAPPER build (osx)".green
    Rake::Task[:build_osx].invoke
  when :linux
    puts "LAPACK_WRAPPER build (linux)".green
    Rake::Task[:build_linux].invoke
  when :win
    puts "LAPACK_WRAPPER build (windows)".green
    Rake::Task[:build_win].invoke
  when :mingw
    puts "LAPACK_WRAPPER build (mingw)".green
    Rake::Task[:build_mingw].invoke
  end
end

desc "clean LAPACK_WRAPPER"
task :clean do
  case OS
  when :mac
    puts "LAPACK_WRAPPER clean (osx)".green
    Rake::Task[:clean_osx].invoke
  when :linux
    puts "LAPACK_WRAPPER clean (linux)".green
    Rake::Task[:clean_linux].invoke
  when :win
    puts "LAPACK_WRAPPER clean (windows)".green
    Rake::Task[:clean_win].invoke
  when :mingw
    puts "LAPACK_WRAPPER clean (mingw)".green
    Rake::Task[:clean_mingw].invoke
  end
end

task :mkl, [:year, :bits] do |t, args|
  args.with_defaults(:year => "2017", :bits => "x64" )
  sh "'C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/bin/compilervars.bat' -arch #{args.bits} vs#{args.year}shell"
end

#desc "build lib"
#task :build do
#  sh "make config"
#  sh "make --jobs=8 install_local"
#end

desc "run tests"
task :test do
  FileUtils.cd "build"
  sh 'ctest --output-on-failure'
  FileUtils.cd '..'
end

TESTS = [
  "test1-small-factorization",
  "test2-Timing",
  "test3-BandedMatrix",
  "test4-BFGS",
  "test5-BLOCKTRID",
  "test6-EIGS",
  "test7-SparseTool",
  "test8-SparseToolComplex",
  "test9-SparseToolTridiagonal",
  "test10-SparseToolTiming",
  "test11-SparseToolVector",
  "test12-SparseToolMatrix",
  "test13-SparseTool1",
  "test14-SparseTool2"
];

desc "run tests"
task :run do
  puts "UTILS run tests".green
  case OS
  when :mac,:linux
    TESTS.each do |cmd|
      exe = "./bin/#{cmd}"
      next unless File.exist?(exe)
      puts "execute #{exe}".yellow
      sh exe
    end
  when :win
    TESTS.each do |cmd|
      exe = "bin\\#{cmd}.exe"
      next unless File.exist?(exe)
      puts "execute #{exe}".yellow
      sh exe
    end
  end
end
def ChangeOnFile( file, text_to_replace, text_to_put_in_place )
  text = File.read file
  File.open(file, 'w+'){|f| f << text.gsub(text_to_replace, text_to_put_in_place)}
end


desc "compile for Visual Studio [default year=2017, bits=x64, lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :build_win, [:year, :bits, :lapack] do |t, args|
  args.with_defaults(
    :year   => "2017",
    :bits   => "x64",
    :lapack => "LAPACK_WRAPPER_USE_OPENBLAS"
  )
  Rake::Task[:win_3rd].invoke(args.year,args.bits,args.lapack)
  Rake::Task[:build_win_common].invoke(args.year,args.bits,args.lapack)
end


desc "compile for Visual Studio [default year=2017, bits=x64, lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :build_mingw do
  Rake::Task[:mingw_3rd].invoke()
  Rake::Task[:build_win_common].invoke()
end


desc "compile for Visual Studio [default year=2017, bits=x64, lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :build_win_common, [:year, :bits, :lapack] do |t, args|
  args.with_defaults(
    :year   => "2017",
    :bits   => "x64",
    :lapack => "LAPACK_WRAPPER_USE_OPENBLAS"
  )

  cmd = "set path=%path%;lib3rd\\lib;lib3rd\\dll;"

  FileUtils.rm_f 'src/lapack_wrapper_config.hh'
  FileUtils.cp   'src/lapack_wrapper_config.hh.tmpl', 'src/lapack_wrapper_config.hh'

  ChangeOnFile(
    'src/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_USE@@',
    "#define #{args.lapack} 1"
  )
  ChangeOnFile(
    'src/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_NOSYSTEM_OPENBLAS@@',
    "#define LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS 1"
  )

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  FileUtils.mkdir_p "../lib/lib"
  FileUtils.mkdir_p "../lib/bin"
  FileUtils.mkdir_p "../lib/bin/"+args.bits
  FileUtils.mkdir_p "../lib/dll"
  FileUtils.mkdir_p "../lib/include"

  cmd_cmake = win_vs(args.bits,args.year) + cmd_cmake_build + ' -D' + args.lapack + ':VAR=ON '

  puts "run CMAKE for LAPACK WRAPPER".yellow

  sh cmd_cmake + ' ..'
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end

desc "compile for Visual Studio [default year=2017, bits=x64, lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :build_mingw do

end

desc 'compile for OSX [default lapack="LAPACK_WRAPPER_USE_ACCELERATE"]'
task :build_osx, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_ACCELERATE" )

  FileUtils.rm_f 'src/lapack_wrapper_config.hh'
  FileUtils.cp   'src/lapack_wrapper_config.hh.tmpl', 'src/lapack_wrapper_config.hh'

  Rake::Task[:osx_3rd].invoke()

  ChangeOnFile(
    'src/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_USE@@',
    "#define #{args.lapack} 1"
  )
  ChangeOnFile(
    'src/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_NOSYSTEM_OPENBLAS@@',
    "// #define LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS 1"
  )

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = "cmake " + cmd_cmake_build + ' -D' + args.lapack + ':VAR=ON '

  puts "run CMAKE for LAPACK WRAPPER".yellow

  sh cmd_cmake + ' ..'
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  if RUN_CPACK then
    puts "run CPACK for SPLINES".yellow
    sh 'cpack -C CPackConfig.cmake'
    sh 'cpack -C CPackSourceConfig.cmake'
  end

  FileUtils.cd '..'
end

desc 'compile for LINUX [default lapack="LAPACK_WRAPPER_USE_OPENBLAS"]'
task :build_linux, [:lapack] do |t, args|
  args.with_defaults(
    :lapack => "LAPACK_WRAPPER_USE_OPENBLAS"
    #:lapack => "LAPACK_WRAPPER_USE_LAPACK",
    #:lapack => "LAPACK_WRAPPER_USE_MKL"
  )

  Rake::Task[:linux_3rd].invoke()

  FileUtils.rm_f 'src/lapack_wrapper_config.hh'
  FileUtils.cp   'src/lapack_wrapper_config.hh.tmpl', 'src/lapack_wrapper_config.hh'

  ChangeOnFile(
    'src/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_USE@@',
    "#define #{args.lapack} 1"
  )
  ChangeOnFile(
    'src/lapack_wrapper_config.hh',
    '@@LAPACK_WRAPPER_NOSYSTEM_OPENBLAS@@',
    "// #define LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS 1"
  )

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = "cmake " + cmd_cmake_build + ' -D' + args.lapack + ':VAR=ON '

  puts "run CMAKE for LAPACK WRAPPER".yellow

  sh cmd_cmake + ' ..'
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  if RUN_CPACK then
    puts "run CPACK for SPLINES".yellow
    sh 'cpack -C CPackConfig.cmake'
    sh 'cpack -C CPackSourceConfig.cmake'
  end

  FileUtils.cd '..'
end

desc 'install third parties for osx'
task :osx_3rd do
  FileUtils.cd 'ThirdParties'
  sh "rake install_osx"
  FileUtils.cd '..'
end

desc 'install third parties for linux'
task :linux_3rd do
  FileUtils.cd 'ThirdParties'
  sh "rake install_linux"
  FileUtils.cd '..'
end

desc "compile for Visual Studio [default year=2017, bits=x64, lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :win_3rd, [:year, :bits, :lapack] do |t, args|
  args.with_defaults(
    :year   => "2017",
    :bits   => "x64",
    :lapack => "LAPACK_WRAPPER_USE_OPENBLAS"
  )
  FileUtils.cd 'ThirdParties'
  sh "rake install_win[#{args.year},#{args.bits},#{args.lapack}]"
  FileUtils.cd '..'
end

desc "compile for Visual Studio [default year=2017, bits=x64, lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :mingw_3rd do
  FileUtils.cd 'ThirdParties'
  sh "rake install_mingw"
  FileUtils.cd '..'
end

desc "clean for OSX"
task :clean_osx do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  FileUtils.rm_rf 'build'
  FileUtils.cd 'ThirdParties'
  sh "rake clean_osx"
  FileUtils.cd '..'
end

desc "clean for LINUX"
task :clean_linux do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  FileUtils.rm_rf 'build'
  FileUtils.cd 'ThirdParties'
  sh "rake clean_linux"
  FileUtils.cd '..'
end

desc "clean for WINDOWS"
task :clean_win do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  FileUtils.rm_rf 'vs_*'
  FileUtils.cd 'ThirdParties'
  sh "rake clean_win"
  FileUtils.cd '..'
end

desc "clean for MINGW"
task :clean_win do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  FileUtils.rm_rf 'build'
  FileUtils.cd 'ThirdParties'
  sh "rake clean_mingw"
  FileUtils.cd '..'
end

task :cppcheck do
  FileUtils.rm_rf   'lib'
  FileUtils.rm_rf   'build'
  FileUtils.mkdir_p 'build'
  FileUtils.cd      'build'
  sh 'cmake -DCMAKE_EXPORT_COMPILE_COMMAND=ON ..'
  sh 'cppcheck --project=compile_commands.json'
end

