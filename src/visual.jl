function visualize_sampling(input::Input, output)

    Drawing(1500, 1500, "basic_test.png")
    origin()
    background("white")
    setopacity(1)

    setline(8)
    sethue("black")
    polysmooth(box(O, 1000, 720, vertices=true), 10, :stroke)
    fontsize(120)
    fontface("italic")
    textcentred("U")

    setline(5)
    inter = 700/(input.m-1)
    counter = 1

    for i in 10:inter:710
        line(Point(-550,i-360), Point(-500,i-360), :stroke)
        line(Point(500,i-360), Point(550,i-360), :stroke)

        if input.r.state[counter] == 1
            p = Point(-550,i-360)
            sethue("blue")
            circle(p,15,:fill)
            sethue("black")
        end

        counter += 1
    end

    sethue("red")
    fontsize(25)
    fontface("italic")
    counter = countmap(output)
    for i = 1:length(output)
        p = Point(550, output[i]*inter-350-inter)
        circle(p,15,:fill)
        if counter[output[i]] > 1
            sethue("black")
            Luxor.text(string(counter[output[i]]), p, valign=:center, halign=:center)
            sethue("red")
        end
    end

    clipreset()
    finish()
    preview()
end

function visualize_proba(input::Input, output, data)

    nlist = output_mode_occupation(input.n, input.m)

    if output in nlist
        idx = findfirst(x -> x==output, nlist)
        print("event probability: ", data[idx])
        visualize_sampling(input, output)
    else
        throw(ArgumentError("invalid argument(s)"))
    end

end
