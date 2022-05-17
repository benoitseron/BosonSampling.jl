# search for counter examples to bapat sunder to show that it is hard
"""
    search_until_user_stop(search_function)

Runs `search_function` until user-stop (Ctrl+C).
"""
function search_until_user_stop(search_function)

    n_trials = 1

    try

        while true
            search_function()
            n_trials += 1
        end

    catch err

        if isa(err, InterruptException)
            print("no counter example found after ")
        end
    finally
        println("$n_trials trials")
    end
end

"""
    random_search_counter_example_bapat_sunder(;m,n,r, physical_H = true)

Brute-force search of counter-examples of rank `r`.
"""
function random_search_counter_example_bapat_sunder(;m,n,r, physical_H = true)

    S = rand_gram_matrix_rank(n,r)

    if physical_H
        U = rand_haar(m)

        input_state = [i <= n ? 1 : 0 for i in 1:m]
        partition = [i <= r ? 1 : 0 for i in 1:m]

        H = H_matrix(U, input_state, partition)
    else
        H = rand_gram_matrix_rank(n)
    end

    if violates_bapat_sunder(H,S)
        println("counter example found : ")
        println(H)
        println(S)
    end


end
