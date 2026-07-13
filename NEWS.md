# nat.ggplot 1.1.6

* `ggneuron_gif()` gains a `points` argument (with `point_cols`, `point_size`,
  `point_alpha`, `point_stroke`) to draw point sets — e.g. neuron root points /
  somata — as circles on top of the neurons. Each entry may be a static matrix or,
  like `flows`, a per-timepoint list so the circles move with the warp.
