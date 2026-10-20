ns_fn <- function(name) {
  getFromNamespace(name, "MPCurver")
}

cavi <- ns_fn("cavi")
do_cavi <- ns_fn("do_cavi")
soft_partition_cavi <- ns_fn("soft_partition_cavi")
as_mpcurve <- ns_fn("as_mpcurve")
init_m_trajectories_cavi <- ns_fn("init_m_trajectories_cavi")
