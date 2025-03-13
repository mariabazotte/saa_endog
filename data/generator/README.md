# Instance File

This README provides details on the structure of the instance files and the parameter values used to generate them. The instance files referenced in this document are contained in the `data/usa` folder and are used in the associated paper.

## Organization of the Instance Files

### First Line
The first line of each instance file contains five integer values representing:
- The number of nodes
- The number of edges
- The number of facilities
- The number of disruption types
- The number of disruption levels

Format:
```
nb_nodes nb_edges nb_facilities nb_disrup_types nb_disrup_levels
```

### Node Information (One Line per Node)
Each subsequent line represents a node, providing the following details:
- The node name (i.e., the city it represents)
- The node’s demand
- The disruption level at which this node is impacted for each disruption type (a value between `0` and `nb_disrup_levels-1`)
- A binary indicator specifying whether the node is a facility (`1` for a facility, `0` otherwise)

Format:
```
name_node demand_node disrup_level_for_disrupt_type_1 ... disrup_level_for_disrupt_type_nb_disrup_types node_is_facility
```

### Edge Information (One Line per Edge)
Each line describing an edge includes:
- The names of the two nodes it connects
- The edge capacity

Format:
```
name_node_source name_node_dest edge_capacity
```

### Disruption Type Information (One Line per Disruption Type)
Each line contains information about a disruption type, including:
- The node that represents the center of the disruption
- The probability of this disruption occurring

Format:
```
name_node_center_disrup_type probab_disrup_type
```

## Parameters Used to Generate Instances in the Paper
The instance files were generated using `usa_network_generator.py`, as described in the main README of this repository. The script requires the number of facilities and a minimum population threshold. Only cities with a population above the specified threshold were selected for inclusion in the network.

The parameters used for different instance sizes are as follows:

| Number of Nodes | Minimum Population |
|----------------|--------------------|
| 15            | 650,000            |
| 25            | 403,500            |
| 30            | 350,000            |
| 39            | 300,000            |
| 48            | 250,000            |

A base script named `job_instance.sh` with these parameters is provided in this folder.