public_exploratory_fixture <- function() {
  data.frame(
    fNode_Coords = c("1_1_1", "2_2_2", "3_3_3", "1_2_1", "2_1_2"),
    tracks = c(4L, 2L, 0L, 1L, 1L),
    local = c(1L, 0L, 0L, 1L, 0L),
    leave_no_target = c(0L, 0L, 0L, 0L, 1L),
    A = c(2L, 2L, 0L, 0L, 0L),
    B = c(2L, 0L, 0L, 0L, 0L), none = c(0L, 0L, 1L, 0L, 0L)
  )
}
