function runs_without_errors(f::Function)

    @test begin
        try
            f()
            true
        catch
            false
        end
    end

end
