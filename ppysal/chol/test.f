C File test.f
    subroutine f1(n)
        integer n
        do 100 i=0, n
            print *, 'Working'
100 continue
    end