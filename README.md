phylogenetic_tree_parser
========================
This project is for my undergraduate thesis. It is an algorithm design and simulation study focusing on phylogentic tree parser. I designed a new phylogenetic tree parsing algorithm from statistical approach with a simple substitution model, JC69 (Jukes and Cantor, 1969) and compared it with a wild used phylogenetic tree parser, Neighbor-Joining (Saitou and Nei, 1987), method with a simulation study on a designed rooted tree with 8 species. 

                        |--------1------sp4
                        |
    |----------1--------|               |---0.01---sp5
    |                   |               |
    |                   |--------1------|           |---0.01----sp6
    |                                   |           |
    |          |-------1-------sp3      |---0.01----| 
    |          |                                    |            |---------1------sp7
    |---0.01---|                                    |---0.01-----|
               |                                                 |---------1------sp8
               |          |---0.01--sp2
               |---0.01---|
                          |----------1---------sp1

