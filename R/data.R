#' LHPD2a1 Neuron Skeletons from BANC Connectome
#'
#' @description
#' A neuronlist containing skeletonised morphologies of 5 LHPD2a1 neurons
#' from the BANC (Brain And Nerve Cord) connectome dataset.
#' These neurons are lateral horn local neurons first described as LHPD2a1 by
#' Dolan et al. (2018) in the adult fly brain. The skeletons have been
#' processed using wavefront skeletonisation and rerooted.
#'
#' @format A neuronlist object containing 5 neurons with 3D point coordinates
#' and connectivity information.
#'
#' @source BANC connectome dataset (Bates et al., 2025). Data accessed via the
#' \code{bancr} package.
#'
#' @references
#' Dolan, M. J., Belliart-Guérin, G., Bates, A. S., Frechter, S., Lampin-Saint-Amaux, A.,
#' Aso, Y., Roberts, R. J. V., Schlegel, P., Wong, A., Hammad, A., Bock, D.,
#' Rubin, G. M., Preat, T., Plaçais, P. Y., & Jefferis, G. S. X. E. (2018).
#' Communication from learned to innate olfactory processing centers is required
#' for memory retrieval in Drosophila. \emph{Neuron}, 100(3), 651-668.e8.
#' \doi{10.1016/j.neuron.2018.08.037}
#'
#' Bates, A. S., Phelps, J. S., Kim, M., Yang, H. H., Matsliah, A., Ajabi, Z.,
#' Perlman, E., Dunne, J., Roat, J., Joyce, P., Bogovic, J. A., Jefferis, G. S. X. E.,
#' Murthy, M., Card, G., & The FlyWire Consortium. (2025).
#' Distributed control circuits across a brain-and-cord connectome. \emph{bioRxiv}.
#' \doi{10.1101/2025.07.31.667571}
#'
#' @examples
#' \dontrun{
#' library(nat)
#' library(nat.ggplot)
#'
#' # Plot the neurons
#' plot3d(banc.skels)
#'
#' # Use with ggplot2
#' library(ggplot2)
#' ggplot() + geom_neuron(banc.skels, rotation_matrix = banc_view)
#' }
"banc.skels"

#' LHPD2a1 Neuron Synapses from BANC Connectome
#'
#' @description
#' A data frame containing synapse locations and metadata for 5 LHPD2a1 neurons
#' from the BANC connectome dataset. Each row represents a single synapse with
#' its 3D position and pre/post synaptic classification.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{connector_id}{Unique identifier for the synapse connector}
#'   \item{pre_id}{Pre-synaptic neuron ID}
#'   \item{post_id}{Post-synaptic neuron ID}
#'   \item{prepost}{Binary indicator: 0 for presynaptic (output), 1 for postsynaptic (input)}
#'   \item{pre_svid}{Pre-synaptic supervoxel ID}
#'   \item{post_svid}{Post-synaptic supervoxel ID}
#'   \item{size}{Size of the synapse}
#'   \item{X, Y, Z}{3D coordinates of the synapse in nanometers}
#'   \item{treenode_id}{Node ID on the skeleton where this synapse is located}
#'   \item{root_id}{The root ID of the neuron this synapse belongs to}
#' }
#'
#' @source BANC connectome dataset (Bates et al., 2025). Data accessed via the
#' \code{bancr} package.
#'
#' @references
#' Bates, A. S., Phelps, J. S., Kim, M., Yang, H. H., Matsliah, A., Ajabi, Z.,
#' Perlman, E., Dunne, J., Roat, J., Joyce, P., Bogovic, J. A., Jefferis, G. S. X. E.,
#' Murthy, M., Card, G., & The FlyWire Consortium. (2025).
#' Distributed control circuits across a brain-and-cord connectome. \emph{bioRxiv}.
#' \doi{10.1101/2025.07.31.667571}
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(nat.ggplot)
#'
#' # Plot synapses coloured by type
#' ggplot(banc.syns, aes(x = X, y = Y)) +
#'   geom_point(aes(colour = factor(prepost)), alpha = 0.5) +
#'   scale_colour_manual(values = c("0" = "#D72000", "1" = "#132157"),
#'                      labels = c("Presynaptic", "Postsynaptic")) +
#'   coord_fixed()
#' }
"banc.syns"

#' Flow Centrality Split LHPD2a1 Neurons from BANC
#'
#' @description
#' A neuronlist containing LHPD2a1 neurons that have been split into axonal
#' and dendritic compartments using flow centrality analysis. This method,
#' described by Schneider-Mizell et al. (2016), uses the spatial distribution
#' of pre- and postsynaptic sites to segregate neuronal arbors into functionally
#' distinct compartments.
#'
#' @format A neuronlist object with neurons containing additional Label field in
#' their data frames indicating compartment type:
#' \describe{
#'   \item{Label = 2}{Axon}
#'   \item{Label = 3}{Dendrite}
#'   \item{Label = 4}{Primary dendrite}
#'   \item{Label = 7}{Primary neurite}
#'   \item{Label = 0 or NA}{Unclassified}
#' }
#'
#' @source BANC connectome dataset (Bates et al., 2025). Split using flow centrality
#' method from \code{hemibrainr} package.
#'
#' @references
#' Schneider-Mizell, C. M., Gerhard, S., Longair, M., Kazimiers, T., Li, F.,
#' Zwart, M. F., Champion, A., Midgley, F. M., Fetter, R. D., Saalfeld, S.,
#' & Cardona, A. (2016). Quantitative neuroanatomy for connectomics in Drosophila.
#' \emph{eLife}, 5, e12059. \doi{10.7554/eLife.12059}
#'
#' Bates, A. S., Phelps, J. S., Kim, M., Yang, H. H., Matsliah, A., Ajabi, Z.,
#' Perlman, E., Dunne, J., Roat, J., Joyce, P., Bogovic, J. A., Jefferis, G. S. X. E.,
#' Murthy, M., Card, G., & The FlyWire Consortium. (2025).
#' Distributed control circuits across a brain-and-cord connectome. \emph{bioRxiv}.
#' \doi{10.1101/2025.07.31.667571}
#'
#' @examples
#' \dontrun{
#' library(nat.ggplot)
#' library(ggplot2)
#'
#' # Plot split neurons showing axon/dendrite compartments
#' ggplot() + geom_neuron(banc.neurons.flow[[1]], rotation_matrix = banc_view)
#'
#' # Plot all split neurons
#' ggneuron(banc.neurons.flow, rotation_matrix = banc_view)
#' }
"banc.neurons.flow"

#' Low Resolution BANC Brain Neuropil Mesh
#'
#' @description
#' A simplified mesh3d object representing the brain neuropil volume from the
#' BANC connectome dataset. This mesh has been decimated to approximately 50%
#' of its original face count for efficient visualisation.
#'
#' @format A mesh3d object containing vertices and faces defining the 3D surface
#' of the Drosophila brain neuropil.
#'
#' @source BANC connectome dataset, decimated using \code{Rvcg::vcgQEdecim}.
#'
#' @references
#' Bates, A. S., Phelps, J. S., Kim, M., Yang, H. H., Matsliah, A., Ajabi, Z.,
#' Perlman, E., Dunne, J., Roat, J., Joyce, P., Bogovic, J. A., Jefferis, G. S. X. E.,
#' Murthy, M., Card, G., & The FlyWire Consortium. (2025).
#' Distributed control circuits across a brain-and-cord connectome. \emph{bioRxiv}.
#' \doi{10.1101/2025.07.31.667571}
#'
#' @examples
#' \dontrun{
#' library(nat.ggplot)
#' library(ggplot2)
#'
#' # Plot brain neuropil with default view
#' ggplot() +
#'   geom_neuron(banc.brain_neuropil,
#'               cols = c("grey75", "grey50"),
#'               alpha = 0.3)
#'
#' # Plot with frontal view
#' ggplot() +
#'   geom_neuron(banc.brain_neuropil,
#'               rotation_matrix = banc_view,
#'               cols = c("grey75", "grey50"),
#'               alpha = 0.3)
#' }
"banc.brain_neuropil"

#' Low Resolution LHPD2a1 Neuron Meshes from BANC
#'
#' @description
#' A list of simplified mesh3d objects representing the 3D surface reconstructions
#' of 4 LHPD2a1 neurons from the BANC connectome dataset.
#'
#' @format A list of 5 mesh3d objects, each containing vertices and faces defining
#' the 3D surface of a neuron.
#'
#' @source BANC connectome dataset, accessed via \code{bancr::banc_read_neuron_meshes}
#' and decimated using \code{Rvcg::vcgQEdecim}.
#'
#' @references
#' Dolan, M. J., Belliart-Guérin, G., Bates, A. S., Frechter, S., Lampin-Saint-Amaux, A.,
#' Aso, Y., Roberts, R. J. V., Schlegel, P., Wong, A., Hammad, A., Bock, D.,
#' Rubin, G. M., Preat, T., Plaçais, P. Y., & Jefferis, G. S. X. E. (2018).
#' Communication from learned to innate olfactory processing centers is required
#' for memory retrieval in Drosophila. \emph{Neuron}, 100(3), 651-668.e8.
#' \doi{10.1016/j.neuron.2018.08.037}
#'
#' Bates, A. S., Phelps, J. S., Kim, M., Yang, H. H., Matsliah, A., Ajabi, Z.,
#' Perlman, E., Dunne, J., Roat, J., Joyce, P., Bogovic, J. A., Jefferis, G. S. X. E.,
#' Murthy, M., Card, G., & The FlyWire Consortium. (2025).
#' Distributed control circuits across a brain-and-cord connectome. \emph{bioRxiv}.
#' \doi{10.1101/2025.07.31.667571}
#'
#' @examples
#' \dontrun{
#' library(nat.ggplot)
#' library(ggplot2)
#'
#' # Plot first neuron mesh
#' ggplot() +
#'   geom_neuron(banc.meshes[[1]],
#'               rotation_matrix = banc_view,
#'               cols = c("purple", "magenta"))
#'
#' # Plot all neuron meshes
#' ggneuron(banc.meshes,
#'          rotation_matrix = banc_view,
#'          cols = c("purple", "magenta"))
#' }
"banc.meshes"

#' BANC Frontal View Rotation Matrix
#'
#' @description
#' A 4x4 rotation matrix for displaying BANC connectome data in a frontal view.
#' This matrix can be applied to any 3D neuroanatomy data to rotate it into a
#' standard frontal viewing orientation. Provded as an example for the user.
#'
#' @format A 4x4 numeric matrix representing a 3D rotation transformation.
#'
#' @source Derived from \code{bancr:::banc_rotation_matrices[["front"]]}.
#'
#' @details
#' Users can create their own custom rotation matrices interactively using the
#' \code{rgl_view()} function after positioning neurons in their desired
#' orientation with \code{plot3d()}.
#'
#' @examples
#' \dontrun{
#' library(nat.ggplot)
#' library(nat)
#' library(ggplot2)
#'
#' # Use the frontal view for plotting
#' ggplot() +
#'   geom_neuron(banc.skels, rotation_matrix = banc_view)
#'
#' # Create a custom view interactively
#' plot3d(banc.brain_neuropil)
#' # Rotate to desired angle with mouse
#' my_view <- rgl_view()$userMatrix
#' # Use custom view
#' ggplot() +
#'   geom_neuron(banc.brain_neuropil, rotation_matrix = my_view)
#' }
"banc_view"
