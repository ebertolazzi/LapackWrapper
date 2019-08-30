#!/usr/bin/env ruby -w

require 'fileUtils'

$figname  = ""
$figlines = ""

def write_figure()
  # stampo comando di inclusione figura
  dir      = File.dirname(__FILE__)+"/DoxygenImages/" ;
  basename = dir+$figname.to_s;
  fname    = basename+".asy";
  outname  = basename+".pdf";
  docmp    = File.exist?(fname) # controllo se il file esiste gia
  File.rename(fname, fname+"_tmp") if docmp
  File.open(fname,'w') do |outfile|
    outfile.puts $figlines
  end
  # controlla se è cambiato il file
  if docmp then
    docmp = FileUtils.compare_file(fname+"_tmp", fname)
    File.delete(fname+"_tmp")
  end
  # docmp = true se il file c'è gia ed è identico al precedente
  if not docmp then
    # il file è cambiato ricompilo
    if !system( "asy -tex pdflatex #{fname} -o #{outname}" ) then
      $stderr.puts "\n\n\n\nasy for #{fname} failed!!!!\n\n\n\n\n"
    end
    system( "pdf2ps #{basename}.pdf #{basename}.eps" ) ;
    system( "convert -density 200 -units PixelsPerInch #{basename}.pdf #{basename}.png" ) ;
  end
  STDOUT.puts "  \\image html #{$figname}.png\n"
  STDOUT.puts "  \\image latex #{$figname}.eps width=10cm\n"
end

$state = :none
#$stderr.puts ARGV[0]
#$stderr.puts 'LANG', ENV['LANG']
#$stderr.puts 'LC_ALL', ENV['LC_ALL']
File.open(ARGV[0],"r").each_line do |line|
  case $state
  when :none
    if line =~ /^\s*(\*\s*)?\\begin\s*\{asy\}\s*\{(\w+)\}/ then
      puts
      $figname  = $2
      $figlines = "locale(\"en_US\");\n"
      $state = :infigure
    else
      STDOUT.puts line
    end
  when :infigure
    if line =~ /^\s*(\*\s*)?\\end\s*\{asy\}/ then
      $state = :none
      write_figure()
    else
      if line =~ /^\s*\*\s?(.*)$/ then
        $figlines += $1 + "\n" ;
      else
        $figlines += line
      end
    end
  end
end

