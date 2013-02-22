context('Correct input and output dimensions and dimension names')

test_that('only original; only K_1; no L; default radius vector', {
    pattern <- spatstat::spruces

    res_l <- summ_func(pattern, mtf_name = '1', do_besags_L = FALSE)
    a <- res_l[['a']]
    r <- res_l[['r']]

    expect_that(a, is_a('matrix'))
    expect_that(dim(a), is_identical_to(c(1L, length(r))))
    expect_that(names(dimnames(a)), is_identical_to(c('summ_func', 'r')))
})

test_that('only original; only K_1; no L; one radius value', {
    pattern <- spatstat::spruces
    r_vec <- 5

    res_l <- summ_func(pattern, mtf_name = '1', r_vec = r_vec,
                       do_besags_L = FALSE)
    a <- res_l[['a']]

    expect_that(a, is_a('matrix'))
    expect_that(dim(a), is_identical_to(c(1L, length(r_vec))))
    expect_that(names(dimnames(a)), is_identical_to(c('summ_func', 'r')))
})

test_that('only original; only K_1; include L; default radius vector', {
    pattern <- spatstat::spruces

    res_l <- summ_func(pattern, mtf_name = '1', do_besags_L = TRUE)
    a <- res_l[['a']]
    r <- res_l[['r']]

    expect_that(a, is_a('matrix'))
    expect_that(dim(a), is_identical_to(c(2L, length(r))))
    expect_that(names(dimnames(a)), is_identical_to(c('summ_func', 'r')))
})



test_random_labelling_dim <- function(pattern, n_perm, mtf_name,
                                      do_besags_L, r_vec) {
    res_l <- summ_func_random_labelling(pattern, n_perm = n_perm,
                                        mtf_name = mtf_name,
                                        do_besags_L = do_besags_L,
                                        r_vec = r_vec)
    a <- res_l[['a']]
    r <- res_l[['r']]

    if (do_besags_L) {
        marks <- pattern[['marks']]
        n_L <- length(besags_L_valid(mtf_name,
                                     is_any_mark_neg = any(marks < 0),
                                     is_any_mark_pos = any(marks > 0)))
    } else {
        n_L <- 0L
    }

    if (length(r_vec) < 1L) {
        n_r <- length(r)
    } else {
        n_r <- length(r_vec)
    }

    expect_that(a, is_a('array'))
    expect_that(dim(a), is_identical_to(c(n_perm + 1L,
                                          length(mtf_name) + n_L,
                                          n_r)))
    expect_that(names(dimnames(a)),
                is_identical_to(c('orig_and_perm', 'summ_func', 'r')))
}

test_that('zero permutations; only K_1; no L; default radius vector', {
    pattern <- spatstat::spruces
    n_perm <- 0L
    mtf_name <- '1'
    do_besags_L <- FALSE
    r_vec <- NULL

    test_random_labelling_dim(pattern, n_perm, mtf_name, do_besags_L, r_vec)
})

test_that('zero permutations; only K_1; no L; one radius value', {
    pattern <- spatstat::spruces
    n_perm <- 0L
    mtf_name <- '1'
    do_besags_L <- FALSE
    r_vec <- 5

    test_random_labelling_dim(pattern, n_perm, mtf_name, do_besags_L, r_vec)
})

test_that('one permutation; only K_1; no L; default radius vector', {
    pattern <- spatstat::spruces
    n_perm <- 1L
    mtf_name <- '1'
    do_besags_L <- FALSE
    r_vec <- NULL

    test_random_labelling_dim(pattern, n_perm, mtf_name, do_besags_L, r_vec)
})

test_that('one permutation; only K_1; no L; one radius value', {
    pattern <- spatstat::spruces
    n_perm <- 1L
    mtf_name <- '1'
    do_besags_L <- FALSE
    r_vec <- 5

    test_random_labelling_dim(pattern, n_perm, mtf_name, do_besags_L, r_vec)
})

test_that('three permutations; only K_1; no L; default radius vector', {
    pattern <- spatstat::spruces
    n_perm <- 3L
    mtf_name <- '1'
    do_besags_L <- FALSE
    r_vec <- NULL

    test_random_labelling_dim(pattern, n_perm, mtf_name, do_besags_L, r_vec)
})

test_that('three permutations; only K_1; no L; one radius value', {
    pattern <- spatstat::spruces
    n_perm <- 3L
    mtf_name <- '1'
    do_besags_L <- FALSE
    r_vec <- 5

    test_random_labelling_dim(pattern, n_perm, mtf_name, do_besags_L, r_vec)
})
