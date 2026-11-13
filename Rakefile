%w(colorize fileutils rake/clean).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end
require 'shellwords'

# Maximum number of parallel compile jobs. GCC/MSVC can need a couple of GB per
# C++20 translation unit, and an over-subscribed build is a classic source of
# compiler crashes. PINS_MAX_CORE (or JOBS) caps it; 16 by default.
def pins_jobs
  require 'etc'
  cap = ENV['PINS_MAX_CORE'].to_i
  cap = 16 if cap <= 0
  n = ENV['JOBS'].to_i
  n = Etc.nprocessors if n <= 0
  [n, cap].min
end


CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }

# Configurazione esterna (opzionale)
if File.exist?(File.expand_path('../Rakefile_configure.rb', File.dirname(__FILE__)))
  require_relative '../Rakefile_configure.rb'
elsif File.exist?(File.expand_path('../../Rakefile_configure.rb', File.dirname(__FILE__)))
  require_relative '../../Rakefile_configure.rb'
else
  COMPILE_DEBUG      = false
  COMPILE_DYNAMIC    = false
  COMPILE_EXECUTABLE = true
end

# ----------------------------------------------------------------------------
# Rilevamento OS
# ----------------------------------------------------------------------------
# MSYSTEM is exported by every MSYS2 shell and is the only reliable way to
# tell a MinGW environment from a plain POSIX one: the ruby shipped in the
# MSYS subsystem reports its platform as "linux"/"cygwin" even when it drives
# a UCRT64 / MINGW64 toolchain, which made this build announce itself as Linux
# and pick the Linux library suffixes. MSYSTEM=MSYS is the POSIX shell and is
# deliberately NOT treated as MinGW.
PINS_MINGW_ENVS = %w(UCRT64 MINGW64 MINGW32 CLANG64 CLANG32 CLANGARM64)
PINS_IS_MINGW   = PINS_MINGW_ENVS.include?(ENV['MSYSTEM'].to_s.upcase)

if PINS_IS_MINGW
  OS = :mingw
else
  case RUBY_PLATFORM
  when /darwin/
    OS = :mac
  when /msys|mingw/
    OS = :mingw
  when /linux|cygwin/
    OS = :linux
  else
    OS = :win
  end
end

if RUBY_PLATFORM =~ /darwin/
  LAPACK_BACKEND = 'Accelerate'
else
  LAPACK_BACKEND = 'OpenBLAS'
end

def build_type
  COMPILE_DEBUG ? 'Debug' : 'Release'
end

def install_prefix
  File.expand_path('lib', __dir__)
end

def project_root
  File.expand_path(__dir__)
end

def build_dir
  File.join(project_root, 'build')
end

def cmake_generator
  'Ninja'
end

def cmake_configure_command(backend, enable_tests:)
  args = [
    "-G", cmake_generator,
    "-B", build_dir,
    "-DCMAKE_BUILD_TYPE=#{build_type}",
    "-DCMAKE_INSTALL_PREFIX=#{install_prefix}",
    "-DBUILD_SHARED_LIBS=#{COMPILE_DYNAMIC ? 'ON' : 'OFF'}",
    "-DBUILD_TESTING=#{enable_tests ? 'ON' : 'OFF'}",
    "-DLAPACK_WRAPPER_INSTALL=ON",
    "-DLAPACK_WRAPPER_STRICT_WARNINGS=OFF",
    "-DLAPACK_WRAPPER_BACKEND=#{backend}",
    project_root
  ]
  args.join(' ')
end

def cmake_build_command(target = nil)
  cmd = ["--build", build_dir, "--config", build_type]
  cmd += ["--target", target] if target
  cmd += ["--parallel", pins_jobs.to_s]
  cmd.join(' ')
end

def configure_cmake(backend, enable_tests:)
  puts "Configuring CMake in #{build_dir} with backend #{backend}...".green
  sh "cmake " + cmake_configure_command(backend, enable_tests: enable_tests)
end

# ----------------------------------------------------------------------------
# Task
# ----------------------------------------------------------------------------
desc "default task --> build"
task :default => :build

desc "git clean reset"
task :git_clean do
  sh "git reset --hard"
  sh "git clean -d -x -f"
end

desc "Configure CMake for library-only build"
task :configure, [:backend] do |t, args|
  backend = args[:backend] || LAPACK_BACKEND
  configure_cmake(backend, enable_tests: false)
end

desc "Configure CMake with tests enabled"
task :configure_tests, [:backend] do |t, args|
  backend = args[:backend] || LAPACK_BACKEND
  configure_cmake(backend, enable_tests: true)
end

desc "Build and install library only"
task :build, [:backend] do |t, args|
  backend = args[:backend] || LAPACK_BACKEND
  configure_cmake(backend, enable_tests: false)
  puts "Compiling and installing library only...".green
  sh "cmake " + cmake_build_command('install')
end

desc "Compile all targets, including tests"
task :compile_all, [:backend] do |t, args|
  backend = args[:backend] || LAPACK_BACKEND
  configure_cmake(backend, enable_tests: true)
  puts "Compiling all targets...".green
  sh "cmake " + cmake_build_command
end

desc "Run tests (configures and compiles tests if needed)"
task :test, [:backend] do |t, args|
  backend = args[:backend] || LAPACK_BACKEND
  configure_cmake(backend, enable_tests: true)
  puts "Compiling tests...".green
  sh "cmake " + cmake_build_command('LapackWrapper_all_tests')
  puts "Running tests...".green
  Dir.chdir(build_dir) do
    sh "ctest -C #{build_type} --output-on-failure"
  end
end

desc "Build and run tests (alias for test)"
task :run, [:backend] => :test do
  # same as test
end

desc "Clean build artifacts (keeps lib/ and lib3rd/)"
task :clean_build do
  puts "Removing build directory...".green
  FileUtils.rm_rf build_dir
end

desc "Full clean (removes build, lib, lib3rd, bin, object files)"
task :clean do
  puts "Full cleaning...".green
  FileUtils.rm_rf build_dir
  FileUtils.rm_rf File.join(project_root, 'lib')
  FileUtils.rm_rf File.join(project_root, 'lib3rd')
  FileUtils.rm_rf File.join(project_root, 'bin')
  FileUtils.rm_f Dir.glob(File.join(project_root, '**', '*.o'))
  FileUtils.rm_f Dir.glob(File.join(project_root, '**', '*.obj'))
end

desc "Package using CPack"
task :cpack do
  puts "Creating packages...".green
  Dir.chdir(build_dir) do
    sh "cpack -C #{build_type} CPackConfig.cmake"
    sh "cpack -C #{build_type} CPackSourceConfig.cmake"
  end
end

# Compatibility aliases.
task build_osx: :build
task build_linux: :build
task build_mingw: :build
task build_win: :build
task clean_osx: :clean
task clean_linux: :clean
task clean_mingw: :clean
task clean_win: :clean

CLOBBER.include []