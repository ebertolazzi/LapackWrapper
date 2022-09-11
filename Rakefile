require_relative "./cmake_utils/Rakefile_common.rb"

CLEAN.include   ["./**/*.o", "./**/*.obj", "./bin/**/example*", "./build"]
CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }
CLOBBER.include []

task :mkl, [:year, :bits] do |t, args|
  args.with_defaults(:year => "2017", :bits => "x64" )
  sh "'C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/bin/compilervars.bat' -arch #{args.bits} vs#{args.year}shell"
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

  FileUtils.mkdir_p "./lib/lib"
  FileUtils.mkdir_p "./lib/bin"
  FileUtils.mkdir_p "./lib/bin/"+args.bits
  FileUtils.mkdir_p "./lib/dll"
  FileUtils.mkdir_p "./lib/include"

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  cmd_cmake = cmake_generation_command(args.bits,args.year) + cmd_cmake_build() + ' -D' + args.lapack + ':VAR=ON '

  puts "run CMAKE for LAPACK WRAPPER".yellow

  sh cmd_cmake + ' ..'
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end

desc "compile for OSX [lapack=LAPACK_WRAPPER_USE_ACCELERATE]"
task :build_osx, [:lapack] do |t, args|
  args.with_defaults(:lapack => "LAPACK_WRAPPER_USE_ACCELERATE")
  Rake::Task[:osx_3rd].invoke()
  Rake::Task[:build_common].invoke(args.lapack)
end

desc "compile for LINUX [lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :build_linux, [:lapack] do |t, args|
  args.with_defaults(:lapack => "LAPACK_WRAPPER_USE_OPENBLAS")
  Rake::Task[:linux_3rd].invoke()
  Rake::Task[:build_common].invoke(args.lapack)
end

desc "compile for MINGW [lapack=LAPACK_WRAPPER_USE_OPENBLAS]"
task :build_mingw, [:lapack] do |t, args|
  args.with_defaults(:lapack => "LAPACK_WRAPPER_USE_OPENBLAS")
  Rake::Task[:mingw_3rd].invoke()
  Rake::Task[:build_common].invoke(args.lapack)
end


desc 'compile for OSX/LINUX/MINGW'
task :build_common, [:lapack] do |t, args|
  args.with_defaults( :lapack => "LAPACK_WRAPPER_USE_ACCELERATE" )

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

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  cmd_cmake = "cmake " + cmd_cmake_build() + ' -D' + args.lapack + ':VAR=ON '

  puts "run CMAKE for LAPACK WRAPPER".yellow

  sh cmd_cmake + ' ..'
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end

desc 'install third parties for OSX'
task :osx_3rd do
  FileUtils.cd 'ThirdParties'
  sh "rake install_osx"
  FileUtils.cd '..'
end

desc 'install third parties for LINUX'
task :linux_3rd do
  FileUtils.cd 'ThirdParties'
  sh "rake install_linux"
  FileUtils.cd '..'
end

desc 'install third parties for MINGW'
task :mingw_3rd do
  FileUtils.cd 'ThirdParties'
  sh "rake install_mingw"
  FileUtils.cd '..'
end

desc 'install third parties for WINDOWS'
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
  FileUtils.rm_rf 'build'
  FileUtils.cd 'ThirdParties'
  sh "rake clean_win"
  FileUtils.cd '..'
end

desc "clean for MINGW"
task :clean_mingw do
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

desc 'pack for OSX/LINUX/WINDOWS'
task :cpack do
  FileUtils.cd "build"
  puts "run CPACK for Embed Figlet".yellow
  sh 'cpack -C CPackConfig.cmake'
  sh 'cpack -C CPackSourceConfig.cmake'
  FileUtils.cd ".."
end
