## unstructured testing!

# this is a sandpit for function testing
# to be converted into proper assurance tests when time permits



#data_surv_rounded_constant <- round_km(data_surv, min_count, method="constant")
#data_surv_rounded_linear <- round_km(data_surv, min_count, method="linear")






approx_safe <- possibly(approx, otherwise=list(y=NA_real_))


data_test <-
  data_surv |>
  # calculate cumulative number of events
  mutate(
    cml.event = cumsum(n.event),
    cml.censor = cumsum(n.censor),
  ) %>%
  # identify groups that have at least `min_count` events
  mutate(

    grp.event1 = pmax(ceiling_any(cml.event, min_count), min_count),
    grp.censor1 = pmax(ceiling_any(cml.censor, min_count), min_count),

    grp.event = cumsum(lag((cml.event%%min_count==0) & n.event>0, 1, 0)),
    grp.censor = cumsum(lag((cml.censor%%min_count==0) & n.censor>0, 1, 0)),

    grp.start.event = ave(x=time, grp.event, FUN = \(x)min(x)-1),
    grp.end.event = ave(x=time, grp.event, FUN = max),
    grp.cml.event = floor(min_count*ave(x=time, grp.event, FUN = cume_dist)) + (ave(x=cml.event, grp.event, FUN = max)-min_count),
    n.event.rounded = grp.cml.event - lag(grp.cml.event,1,0),

    grp.start.censor = ave(x=time, grp.censor, FUN = \(x)min(x)-1),
    grp.end.censor = ave(x=time, grp.censor, FUN = max),
    grp.cml.censor = floor(min_count*ave(x=time, grp.censor, FUN = cume_dist)) + (ave(x=cml.censor, grp.censor, FUN = max)-min_count),
    n.censor.rounded = grp.cml.censor - lag(grp.cml.censor,1,0),

  )









## some test data ----

# example 1
# create some dummy observation times
stime1 <- c(1,1,1,2,2,3,
            4,4,4,4,5,5,
            6,6,6,6,7,7,
            7,7,7,7,7,7,
            7,7,7,8,8,8,
            9,9,9,9,11,11,
            12,13,13)
scount1 <- times2counts(stime1)$n
scmlcount1 <- cumsum(scount1)
stimeindex1 <- times2counts(stime1)$time


# example 2
# create some dummy observation times
scount2 <-   c(0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,1,
               0,0,2,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,
               0,1,1,1,1,0,0,0)
scmlcount2 <- cumsum(scount2)
stime2 <- counts2times(seq_along(scounts2), scounts2)
stimeindex2 <- times2counts(stime2)$time


# run rounding routines on example data
stime <- stime1
scount <- scount1
scmlcount <- scmlcount1
stimeindex <- stimeindex1




## test test

xroundcml <- ceiling_any(scmlcount, min_count)
xroundcount <- diff(c(0,xroundcml))
xcuts <- c(0,unique(stime))[xroundcount>0]
xcounts <- xroundcount[xroundcount>0]

x_floor <- floor_any(scmlcount1, min_count)
x_ceiling <- lag(x_floor, 1, 0L) + min_count
x_ceiling2 <- ceiling_any(scmlcount1, min_count)


# get number of obs
len <- length(stime)
# define minimum count
min_count <- 6
maxgroups <- floor(len/min_count)
maxtime <- max(stime)

# randomly re-order the sample times to check things work when not ordered
stime_scat <- sample(stime)


checkcml1 <- tibble(
  time = stimeindex,
  cml.event = scmlcount,
  n.event = scount,
  mid = roundmid_any(cml.event, min_count),
  ceiling = ceiling_any(cml.event, min_count),
  floor = floor_any(cml.event, min_count),
  group_block0 = lag(floor, 1, 0L),# + min_count,
  group_block1 = lag(floor, 1, 0L) + min_count,

  test_constant = round_cmlcount(cml.event, time, min_count=min_count, method="constant"),
  test_linear = round_cmlcount(cml.event, time, min_count=min_count, method="linear"),
)


check1 <- tibble(
  x = stime,
  x_conjoin_floor = round_event_times(stime, min_count, origin=0, group.method="conjoin", value.method="floor"),
  x_conjoin_ceiling = round_event_times(stime, min_count, origin=0, group.method="conjoin", value.method="ceiling"),
  x_conjoin_mid = round_event_times(stime, min_count, origin=0, group.method="conjoin", value.method="mid"),
  x_conjoin_spaced = round_event_times(stime, min_count, origin=0, group.method="conjoin", value.method="spaced"),
  x_forward_floor = round_event_times(stime, min_count, origin=0, group.method="forward", value.method="floor"),
  x_forward_ceiling = round_event_times(stime, min_count, origin=0, group.method="forward", value.method="ceiling"),
  x_forward_mid = round_event_times(stime, min_count, origin=0, group.method="forward", value.method="mid"),
  x_forward_spaced = round_event_times(stime, min_count, origin=0, group.method="forward", value.method="spaced"),
  x_cumulative_floor = round_event_times(stime, min_count, origin=0, group.method="cumulative", value.method="floor"),
  x_cumulative_ceiling = round_event_times(stime, min_count, origin=0, group.method="cumulative", value.method="ceiling"),
  x_cumulative_mid = round_event_times(stime, min_count, origin=0, group.method="cumulative", value.method="mid"),
  x_cumulative_spaced = round_event_times(stime, min_count, origin=0, group.method="cumulative", value.method="spaced"),
)


checktab1 <-
  check1 |>
  reframe(
    across(
      everything(),
      ~ {
        survival::survfit(survival::Surv(.x, rep(1, length(.x))) ~ 1, conf.type="log-log") |>
          broom::tidy() |>
          select(-std.error, -conf.low, -conf.high, -n.censor, -n.risk, -estimate) |>
          tidyr::complete(
            time = seq_len(maxtime), # fill in 1 row for each day of follow up
            fill = list(n.event = 0L, n.censor = 0L) # fill in zero events on those days
          ) |>
          mutate(
            cml.event = cumsum(n.event)
          )
      }
    )
  )

checkred <- check %>%
  group_by(x) %>%
  mutate(n=n()) %>%
  summarise(across(everything(), first)) %>%
  ungroup() %>%
  mutate(
    cml.n = cumsum(n)
  )

check_long <- check %>%
  pivot_longer(
    cols=everything(),
    names_to="method", values_to="time"
  ) %>%
  arrange(method, time) %>%
  group_by(method) %>%
  summarise(
    index = (0:10),
    ci = (map_dbl(0:10, ~sum(time<=.x))),
  ) %>%
  mutate(
    time_jit = index+runif(n(), -0.2, 0.2),
    ci_jit = ci + runif(n(), -0.2, 0.2),
    group.method = gsub("x\\_(.+)\\_(.+)", "\\1", method),
    value.method = gsub("x\\_(.+)\\_(.+)", "\\2", method),
  )
check_long %>%
  ggplot()+
  geom_step(aes(x=index, y=ci, colour=value.method, group=method))+
  facet_wrap(facets=vars(group.method))+
  theme_bw()


#### flam


test_time_sort <- sort(stime)
rleobj <- rle(test_time_sort)
rlevalues <- rleobj$values
rlelength <- rleobj$lengths
rlecumlength <- cumsum(rlelength)

rlecumceiling <- ceiling_any(rlecumlength, min_count)

counts <- diff(c(0,rlecumceiling))
cuts <- c(0,rlevalues)[counts>0]
counts <- counts[counts>0]



## testing


round_event_times(stime, min_count, origin=0, group.method="conjoin", value.method="floor")

stte <- stime1
sn.event <- times2counts(stte)$n
scml.event <- cumsum(sn.event)

trleobj <- rle(stte)
trletimes <- trleobj$values
trleevents <- trleobj$lengths
trlecmlevents <- cumsum(trleevents)


