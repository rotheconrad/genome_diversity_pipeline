#!/usr/bin/env ruby

eo = ARGV.shift or raise "Usage: #{$0} path/to/eo-file.txt"

def res2hash(txt)
  Hash[txt.strip.split(',').map { |i| i.split('=', 2) }]
end

def size2mb(txt)
  factor = { b: 1e-6, kb: 1e-3, mb: 1, gb: 1e3, tb: 1e6 }
  txt =~ /^(\d+)([tgmk]?b)/ or raise "Cannot parse number: #{txt}"
  $1.to_f * factor[$2.to_sym]
end

def time2min(txt)
  time = txt.split(':').map(&:to_f).reverse
  (time[0] / 60) + (time[1] || 0.0) +
    60.0 * ((time[2] || 0.0) + 24.0 * (time[3] || 0.0))
end

def nodes2cores(txt)
  txt =~ /^(\d+):ppn=(\d+)$/ or raise "Cannot parse nodes: #{txt}"
  $1.to_i * $2.to_i
end

out = [''] * 6
File.open(eo, 'r') do |fh|
  fh.each do |ln|
    if ln =~ /^Resources: +(.+)/
      res = res2hash($1)
      out[0] = size2mb(res['mem'])
      out[2] = time2min(res['walltime'])
      out[4] = out[2] * nodes2cores(res['nodes'])
    elsif ln =~ /^Rsrc Used: +(.+)/
      res = res2hash($1)
      out[1] = size2mb(res['mem'])
      out[3] = time2min(res['walltime'])
      out[5] = time2min(res['cput'])
      break
    end
  end
end if File.exist? eo
puts out.join("\t")

