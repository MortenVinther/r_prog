announce<-function(x) {
  cat("\nobjective:",x$objective,"  convergence:",x$convergence, "  ", x$message, "  iterations:",x$iterations, "  evaluations:",x$evaluations)
}
