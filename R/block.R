#' @export
random_split <- function(train, X_train, k, split.seed = 2019) {
  # Save the old random seed and use the new one, if any
  if (!is.na(split.seed)) {
    if (exists(".Random.seed")) {
      saved.seed <- .Random.seed
    }
    else {
      saved.seed <- NA
    }
    set.seed(split.seed)
  }

  n <- nrow(train)
  nk <- floor(n / k)
  splits <- c(rep(nk, k - 1), n - nk * (k - 1))

  a <- 1:n
  y <- train[, 3]
  coords <- train[, 1:2]

  index.part <- list()
  train.part <- list()
  X_train.part <- list()

  for (i in 1:k) {
    index.part[[i]] <- sample(a, splits[i], replace = FALSE)
    train.part[[i]] <- train[index.part[[i]], ]
    X_train.part[[i]] <- X_train[index.part[[i]], ]
    a <- setdiff(a, index.part[[i]])
  }

  # Restore the old random seed, if any
  if (!is.na(split.seed) && !is.na(saved.seed)) {
    .Random.seed <- saved.seed
  }

  return(list(
    index.part = index.part,
    train.part = train.part,
    X_train.part = X_train.part
  ))
}
