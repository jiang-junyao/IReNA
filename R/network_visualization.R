#' Display the network through igraph package
#'
#' @param network data.frame, indicating network file
#' @param layout the layout to display the network, options: 'grid','sphere','circle','random'
#' @param cols vector, indicating cols of each group
#' @param cex numeric, indicating the font size of node labels
#' @param type network type, 'TF' indicate TFs network, 'module' indicate intramodular network
#' @param legend TRUE or FALSE, whether to show the legend
#' @import ggplot2
#' @importFrom igraph E
#' @importFrom igraph V
#' @importFrom igraph layout_on_grid
#' @importFrom igraph layout_on_sphere
#' @importFrom igraph layout_in_circle
#' @importFrom igraph layout_randomly
#' @export
#' @return figure
#' @examples
#' plot_network(tf_network_cor0.6_thr10, layout='circle', type='TF')
#'
plot_network <- function(network, layout, cols = NULL, cex = NULL, type = "TF", legend = TRUE) {
  if (is.null(cex)) {
    cex1 <- 0.8
  } else {
    cex1 <- cex
  }
  if (is.null(cols)) {
    col1 <- c("OrangeRed", "SlateBlue3", "DarkOrange", "GreenYellow", "Purple", "DarkSlateGray", "Gold", "DarkGreen", "DeepPink2", "Red4", "#4682B4", "#FFDAB9", "#708090", "#836FFF", "#CDC673", "#CD9B1D", "#FF6EB4", "#CDB5CD", "#008B8B", "#43CD80", "#483D8B", "#66CD00", "#CDC673", "#CDAD00", "#CD9B9B", "#FF8247", "#8B7355", "#8B3A62", "#68228B", "#CDB7B5", "#CD853F", "#6B8E23", "#696969", "#7B68EE", "#9F79EE", "#B0C4DE", "#7A378B", "#66CDAA", "#EEE8AA", "#00FF00", "#EEA2AD", "#A0522D", "#000080", "#E9967A", "#00CDCD", "#8B4500", "#DDA0DD", "#EE9572", "#EEE9E9", "#8B1A1A", "#8B8378", "#EE9A49", "#EECFA1", "#8B4726", "#8B8878", "#EEB4B4", "#C1CDCD", "#8B7500", "#0000FF", "#EEEED1", "#4F94CD", "#6E8B3D", "#B0E2FF", "#76EE00", "#A2B5CD", "#548B54", "#BBFFFF", "#B4EEB4", "#00C5CD", "#008B8B", "#7FFFD4", "#8EE5EE", "#43CD80", "#68838B", "#00FF00", "#B9D3EE", "#9ACD32", "#00688B", "#FFEC8B", "#1C86EE", "#CDCD00", "#473C8B", "#FFB90F", "#EED5D2", "#CD5555", "#CDC9A5", "#FFE7BA", "#FFDAB9", "#CD661D", "#CDC5BF", "#FF8C69", "#8A2BE2", "#CD8500", "#B03060", "#FF6347", "#FF7F50", "#CD0000", "#F4A460", "#FFB5C5", "#DAA520", "#CD6889", "#32CD32", "#FF00FF", "#2E8B57", "#CD96CD", "#48D1CC", "#9B30FF", "#1E90FF", "#CDB5CD", "#191970", "#E8E8E8", "#FFDAB9")
  } else {
    col1 <- cols
  }
  if (type == "TF") {
    tfs <- network[, c("TFSymbol", "TFGroup")]
    target <- network[, c("TargetSymbol", "TargetGroup")]
    colnames(target) <- c("TFSymbol", "TFGroup")
    nodes <- rbind(tfs, target)
    edges <- network[, c("TFSymbol", "TargetSymbol", "Regulation", "Correlation")]
  } else if (type == "module") {
    tfs <- network[, c("TFGroup", "TFGroup")]
    target <- network[, c("TargetGroup", "TargetGroup")]
    colnames(target) <- c("TFGroup", "TFGroup")
    nodes <- rbind(tfs, target)
    edges <- network[, c("TFGroup", "TargetGroup", "Regulation", "Correlation")]
  }
  colnames(nodes) <- c("name", "type")
  nodes <- nodes[!duplicated(nodes$name), ]
  colnames(edges) <- c("from", "to", "type", "weight")
  g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
  if (layout == "grid") {
    layout1 <- igraph::layout_on_grid(g)
  } else if (layout == "sphere") {
    layout1 <- igraph::layout_on_sphere(g)
  } else if (layout == "circle") {
    layout1 <- igraph::layout_in_circle(g)
  } else if (layout == "random") {
    layout1 <- igraph::layout_randomly(g)
  } else {
    print("please input correct layout name")
  }
  igraph::E(g)$arrow.size <- 0.3
  igraph::V(g)$color <- col1[nodes$type]
  igraph::V(g)$border <- "red"
  plot(g, layout = layout1, edge.curved = .1, vertex.label.cex = cex1, layout = layout1)
  if (legend == TRUE) {
    legend(x = 1.5, y = 1.3, levels(factor(igraph::V(g)$type)), pch = 21, col = "#777777", pt.bg = col1)
  }
}

#' initiating the cytoscape
#'
#'
#' @param network data.frame, indicating network file
#' @param colour vector, indicating colors for border of each node
#' @param type network type, 'TF' indicate TFs network, 'module' indicate intramodular network
#' @param layout1 the layout to display the network, options: "degree-circle", "attributes-layout", "kamada-kawai", "force-directed", "cose", "hierarchical", "attribute-circle", "stacked-node-layout", "circular", "grid"
#' @return NA
#' @importFrom RCy3 cytoscapePing
#' @importFrom RCy3 cytoscapeVersionInfo
#' @importFrom RCy3 createNetworkFromDataFrames
#' @importFrom RCy3 mapVisualProperty
#' @importFrom RCy3 createVisualStyle
#' @importFrom RCy3 setVisualStyle
#' @importFrom RCy3 layoutNetwork
#' @export
#'
#' @examples \dontrun{
#' initiate_cy(tf_network, layout1='degree-circle', type='TF')
#' initiate_cy(intramodule_network, layout1='grid', type='module')
#' }
initiate_cy <- function(network, colour = NULL, type = "TF", layout1 = "degree-circle") {
  RCy3::cytoscapePing()
  RCy3::cytoscapeVersionInfo()
  if (type == "TF") {
    tfs <- network[, c("TFSymbol", "TFGroup")]
    target <- network[, c("TargetSymbol", "TargetGroup")]
    colnames(target) <- c("TFSymbol", "TFGroup")
    nodes <- rbind(tfs, target)
    edges <- network[, c("TFSymbol", "TargetSymbol", "Regulation", "Correlation")]
  } else if (type == "module") {
    tfs <- network[, c("TFGroup", "TFGroup")]
    target <- network[, c("TargetGroup", "TargetGroup")]
    colnames(target) <- c("TFGroup", "TFGroup")
    nodes <- rbind(tfs, target)
    edges <- network[, c("TFGroup", "TargetGroup", "Regulation", "Correlation")]
    edges$TFGroup <- as.character(edges$TFGroup)
    edges$TargetGroup <- as.character(edges$TargetGroup)
    nodes$TFGroup <- as.character(nodes$TFGroup)
  }
  colnames(nodes) <- c("id", "group")
  nodes <- nodes[!duplicated(nodes$id), ]
  colnames(edges) <- c("source", "target", "interaction", "weight")
  RCy3::createNetworkFromDataFrames(nodes, edges, title = "my first network", collection = "DataFrame Example")
  style.name <- "network1"
  if (is.null(colour)) {
    col1 <- c("OrangeRed", "SlateBlue3", "DarkOrange", "GreenYellow", "Purple", "DarkSlateGray", "Gold", "DarkGreen", "DeepPink2", "Red4", "#4682B4", "#FFDAB9", "#708090", "#836FFF", "#CDC673", "#CD9B1D", "#FF6EB4", "#CDB5CD", "#008B8B", "#43CD80", "#483D8B", "#66CD00", "#CDC673", "#CDAD00", "#CD9B9B", "#FF8247", "#8B7355", "#8B3A62", "#68228B", "#CDB7B5", "#CD853F", "#6B8E23", "#696969", "#7B68EE", "#9F79EE", "#B0C4DE", "#7A378B", "#66CDAA", "#EEE8AA", "#00FF00", "#EEA2AD", "#A0522D", "#000080", "#E9967A", "#00CDCD", "#8B4500", "#DDA0DD", "#EE9572", "#EEE9E9", "#8B1A1A", "#8B8378", "#EE9A49", "#EECFA1", "#8B4726", "#8B8878", "#EEB4B4", "#C1CDCD", "#8B7500", "#0000FF", "#EEEED1", "#4F94CD", "#6E8B3D", "#B0E2FF", "#76EE00", "#A2B5CD", "#548B54", "#BBFFFF", "#B4EEB4", "#00C5CD", "#008B8B", "#7FFFD4", "#8EE5EE", "#43CD80", "#68838B", "#00FF00", "#B9D3EE", "#9ACD32", "#00688B", "#FFEC8B", "#1C86EE", "#CDCD00", "#473C8B", "#FFB90F", "#EED5D2", "#CD5555", "#CDC9A5", "#FFE7BA", "#FFDAB9", "#CD661D", "#CDC5BF", "#FF8C69", "#8A2BE2", "#CD8500", "#B03060", "#FF6347", "#FF7F50", "#CD0000", "#F4A460", "#FFB5C5", "#DAA520", "#CD6889", "#32CD32", "#FF00FF", "#2E8B57", "#CD96CD", "#48D1CC", "#9B30FF", "#1E90FF", "#CDB5CD", "#191970", "#E8E8E8", "#FFDAB9")
  } else {
    col1 <- colour
  }
  defaults <- list(
    NODE_SHAPE = "BioPAX",
    NODE_SIZE = 120,
    EDGE_TRANSPARENCY = 120,
    NODE_LABEL_FONT_SIZE = 29,
    NODE_FILL_COLOR = "White"
  )
  nodeLabels <- RCy3::mapVisualProperty("node label", "id", "p")
  nodeborder <- RCy3::mapVisualProperty(
    "Node Border Paint", "group", "d", levels(as.factor(network$TFGroup)),
    col1
  )
  arrowShapes <- RCy3::mapVisualProperty("Edge Target Arrow Shape", "interaction", "d", levels(as.factor(network$Regulation)), c("Arrow", "T", "None"))
  edgeWidth <- RCy3::mapVisualProperty("edge width", "weight", "p")
  nodeFills <- RCy3::mapVisualProperty("Node Fill Color", "group", "p")
  RCy3::createVisualStyle(style.name, defaults, list(nodeLabels, nodeborder, arrowShapes, edgeWidth, nodeFills))
  RCy3::setVisualStyle(style.name)
  RCy3::layoutNetwork(layout1)
}
