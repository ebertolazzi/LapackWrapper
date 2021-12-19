%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

require_relative "./Rakefile_common.rb"

file_base = File.expand_path(File.dirname(__FILE__)).to_s+'/lib'

cmd_cmake_build = ""
if COMPILE_EXECUTABLE then
  cmd_cmake_build += ' -DBUILD_EXECUTABLE:VAR=true '
else
  cmd_cmake_build += ' -DBUILD_EXECUTABLE:VAR=false '
end
if COMPILE_DYNAMIC then
  cmd_cmake_build += ' -DBUILD_SHARED:VAR=true '
else
  cmd_cmake_build += ' -DBUILD_SHARED:VAR=false '
end
if COMPILE_DEBUG then
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=STATUS '
else
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=STATUS '
end
cmd_cmake_build += " -DEB_INSTALL_LOCAL=ON "


task :default => [:build]

task :mkl, [:year, :bits] do |t, args|
  args.with_defaults(:year => "2017", :bits => "x64" )
  sh "'C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/bin/compilervars.bat' -arch #{args.bits} vs#{args.year}shell"
end

desc "build lib"
task :build do
  sh "make config"
  sh "make --jobs=8 install_local"
end

desc "run tests"
task :test do
  FileUtils.cd "build"
  sh 'ctest --output-on-failure'
  FileUtils.cd '..'
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

  cmd_cmake = win_vs(args.bits,args.year) + cmd_cmake_build

  puts "run CMAKE for LAPACK WRAPPER".yellow

  sh cmd_cmake + ' ..'
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
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

  cmd_cmake = "cmake " + cmd_cmake_build + ' -D' + args.lapack + '=true '

  puts "run CMAKE for LAPACK WRAPPER".yellow

  sh cmd_cmake + ' ..'
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
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

  cmd_cmake = "cmake " + cmd_cmake_build
  # non serve
  # cmd_cmake += ' -D' + args.lapack + '=true '

  puts "run CMAKE for LAPACK WRAPPER".yellow

  sh cmd_cmake + ' ..'
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
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

desc 'pack for OSX/LINUX/WINDOWS'
task :cpack do
  FileUtils.cd "build"
  puts "run CPACK for LAPACK_WRAPPER".yellow
  sh 'cpack -C CPackConfig.cmake'
  sh 'cpack -C CPackSourceConfig.cmake'
  FileUtils.cd ".."
end

desc "clean for OSX"
task :clean_osx do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
  FileUtils.cd 'ThirdParties'
  sh "rake clean_osx"
  FileUtils.cd '..'
end

desc "clean for LINUX"
task :clean_linux do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'lib3rd'
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
