context('spatstat equivalence for r, K, L')

# Test helper functions.

extract_spatstat_result <- function(sp_K_res, sp_L_res, name) {
    list(K = sp_K_res[[name]], L = sp_L_res[[name]])
}

test_own_and_spatstat <- function(pattern, edge_correction, do_besags_L,
                                  r_vec = NULL) {
    our_res <- summ_func(pattern, edge_correction = edge_correction,
                         mtf_name = '1', do_besags_L = do_besags_L,
                         r_vec = r_vec)
    our_r <- our_res[['r']]
    our_a <- our_res[['a']]
    our_K <- our_a['K_1', , drop = TRUE]
    if (do_besags_L) {
        our_L <- our_a['L_1', , drop = TRUE]
    }

    sp_K_res <- spatstat::Kest(pattern, r = r_vec,
                             correction = edge_correction)
    sp_L_res <- NULL
    if (do_besags_L) {
        sp_L_res <- spatstat::Lest(pattern, r = r_vec,
                                   correction = edge_correction)
    }

    if (edge_correction == 'translate') {
        sp_K_L_l <- extract_spatstat_result(sp_K_res, sp_L_res, 'trans')
    } else if (edge_correction == 'none') {
        sp_K_L_l <- extract_spatstat_result(sp_K_res, sp_L_res, 'un')
    }

    sp_K_r <- sp_K_res[['r']]
    sp_L_r <- sp_L_res[['r']]
    sp_K <- sp_K_L_l[['K']]
    sp_L <- sp_K_L_l[['L']]

    expect_that(our_K, equals(sp_K))
    if (length(r_vec) > 0L) {
        expect_that(r_vec, is_identical_to(sp_K_r))
        expect_that(r_vec, is_identical_to(our_r))
    }
    if (do_besags_L) {
        expect_that(sp_K_r, is_identical_to(sp_L_r))
        expect_that(our_L, equals(sp_L))
    }
}



# Tests.

test_that('spruces; r, K and L; default radius vector', {
    pattern <- spatstat::spruces
    test_own_and_spatstat(pattern, edge_correction = 'translate',
                          do_besags_L = TRUE, r_vec = NULL)
})

test_that('spruces; r, K and L; given radius vector', {
    pattern <- spatstat::spruces
    r_vec <- seq(0, 9, length.out = 1e1L)
    test_own_and_spatstat(pattern, edge_correction = 'translate',
                          do_besags_L = TRUE, r_vec = r_vec)
})

test_that('spruces; r, K and L; no edge correction', {
    pattern <- spatstat::spruces
    test_own_and_spatstat(pattern, edge_correction = 'none',
                          do_besags_L = TRUE, r_vec = NULL)
})
