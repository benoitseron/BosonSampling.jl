# Convert julia script files to markdown and then using pandoc to LaTeX
names = filter(name->endswith(name, ".jl") && !occursin("make", name), readdir("."))

# Option specifying whether to execute scripts
run_scripts = length(ARGS) == 0 ? nothing : ARGS[1]

for n=names
  # Get file name
  fn = splitext(n)[1]

  # Run if option is given
  if run_scripts == "run"
    println("Running $fn.jl")
    run(`julia $fn.jl`)
  elseif run_scripts == nothing
    nothing
  else
    error("Unknown option!")
  end


  # Open file
  touch("$fn.md")
  f_in = open(n)
  f_out = open("$fn.md", "w")

  # Add markdown indicator to top
  println(f_out, "```julia")

  # println f_in to f_out
  lines = readlines(f_in)
  for l in lines
     println(f_out, l)
  end
  close(f_in)

  # Add markdown indicator at bottom
  write(f_out, "```")
  close(f_out)

  # Run pandoc for conversion to temporary .tex file
  run(`pandoc $fn.md -o tmp.tex`)

  # Read tmp.tex again and replace the using by keyword highlighting
  touch("$fn.tex")
  tex_in = open("tmp.tex")
  lines2 = readlines(tex_in)
  lines_out = String[]
  hide_block = false
  for l in lines2
      hide_block = hide_block || occursin("hide-block", l)
      if occursin("hide-block-end", l)
          hide_block = false
          continue
      end
      (occursin("hide", l) || hide_block) && continue
      l = replace(l, "\\NormalTok{using" => "\\KeywordTok{using}\\NormalTok{")
      l = replace(l, "\\NormalTok{@cnumbers" => "\\KeywordTok{@cnumbers}\\NormalTok{")
      l = replace(l, "\\NormalTok{@qnumbers" => "\\KeywordTok{@qnumbers}\\NormalTok{")
      l = replace(l, "\\NormalTok{@named" => "\\KeywordTok{@named}\\NormalTok{")
      l = replace(l, "true" => "\\FloatTok{true}")
      l = replace(l, "false" => "\\FloatTok{false}")
      push!(lines_out, l)
  end
  close(tex_in)

  # Delete empty newlines at the end
  while isempty(lines_out[end])
      deleteat!(lines_out, lastindex(lines_out))
  end
  # Delete empty lines before \end{Highlighting}
  while isempty(lines_out[end-2])
      deleteat!(lines_out, lastindex(lines_out)-2)
  end

  tex_out = open("$fn.tex", "w")
  for l in lines_out
      println(tex_out, l)
  end

  close(tex_out)
end

rm("tmp.tex")
