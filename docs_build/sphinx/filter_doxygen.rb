require 'fileutils'

File.open(ARGV[0],"r") do |file|
  file.each_line do |line|
    line.gsub!(/([a-zA-Z]+)\[\]/) { |m| " * "+$1 };
    puts line;
  end
end
