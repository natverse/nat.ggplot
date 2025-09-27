# Create data for package
library(bancr)

# Choose IDs
ids <- c("720575941455119667",
         "720575941554219547",
         "720575941477125441",
         "720575941431479544",
         "720575941483838650")

# Get meta data
banc.meta <- banctable_query()
banc.roots <- banc.meta %>%
  dplyr::distinct(root_id = root_626,
                  root_position,
                  root_position_nm)

# Skeletonise these neurons
banc.skels <- fafbseg::skeletor(segments = ids,
                                cloudvolume.url = bancr:::banc_cloudvolume_url(),
                                clean = TRUE,
                                method = "wavefront",
                                save.obj = NULL,
                                mesh3d = FALSE,
                                waves = 1,
                                k.soma.search = 50,
                                radius.soma.search = 10000,
                                heal = TRUE,
                                reroot = TRUE,
                                heal.threshold = 1000000,
                                heal.k = 10L,
                                reroot_method = "density",
                                brain = bancr::banc_neuropil.surf,
                                elapsed = 100000,
                                resample = 1000)
ids <- names(banc.skels)
meta <- banc.meta %>%
  dplyr::filter(root_626 %in% ids)
rownames(meta) <- meta$root_626

# Get LHPD2a1 neuron meshes
banc.meshes <- bancr::banc_read_neuron_meshes(ids)

# Re-root skeletons
banc.skels <- hemibrainr::add_field_seq(banc.skels, entries=ids, field="id")
banc.skels.rerooted <- nlapply(banc.skels,
                               banc_reroot,
                               id = NULL,
                               roots = banc.roots,
                               estimate = FALSE)

# Add synapses
banc.skels.rerooted <- hemibrainr::add_field_seq(banc.skels.rerooted, entries=ids, field="id")
banc.neurons.syns <- nlapply(banc.skels.rerooted,
                             banc_add_synapses,
                             OmitFailures = TRUE,
                             id = NULL)
banc.neurons.syns <- hemibrainr::add_field_seq(banc.neurons.syns, entries=ids, field="id")
banc.syns <- do.call(rbind, nat::nlapply(banc.neurons.syns, function(x) {
  x$connectors$root_id <- x$id
  x$connectors
}))
rownames(banc.syns) <- 1:nrow(banc.syns)

# Split skeletons
banc.neurons.flow <- hemibrainr::flow_centrality(banc.neurons.syns,
                                               mode = "centrifugal",
                                               polypre = TRUE,
                                               split = "synapses",
                                               .parallel = FALSE,
                                               OmitFailures = FALSE)
banc.neurons.flow[,] <- meta

# # Make a lower res BANC brain mesh
banc.brain_neuropil <- as.mesh3d(bancr::banc_brain_neuropil.surf)
# banc.brain_neuropil <- Rvcg::vcgQEdecim(banc.brain_neuropil, percent = .5)
# # Simplify banc_brain_neuropil and banc.meshes by ~50% face number
# banc.meshes <- nat::nlapply(banc.meshes,  Rvcg::vcgQEdecim, percent = .2)

# # Get BANC nuclei!
# banc.nuclei <- banc_read_nuclei_mesh(meta$nucleus_id)
# banc.nuclei <- nat::nlapply(banc.nuclei,  Rvcg::vcgQEdecim, percent = .2)

# Save meshes: banc.meshes
usethis::use_data(
  banc.meshes,
  overwrite = TRUE,
  compress = "xz"
)

# Save skeletons: banc.neurons.flow
usethis::use_data(
  banc.skels,
  overwrite = TRUE,
  compress = "xz"
)

# Save split skeletons: banc.neurons.flow
usethis::use_data(
  banc.neurons.flow,
  overwrite = TRUE,
  compress = "xz"
)

# Save synapses:
usethis::use_data(
  banc.syns,
  overwrite = TRUE,
  compress = "xz"
)

# Save brain mesh:
usethis::use_data(
  banc.brain_neuropil,
  overwrite = TRUE,
  compress = "xz"
)

# Save brain mesh:
banc_view <- bancr:::banc_rotation_matrices[["front"]]
usethis::use_data(
  banc_view,
  overwrite = TRUE,
  compress = "xz"
)


