/* -*- mode: c++; -*- */
/*
  Prints basic set of statistics about performance of intersection
  accelerators.  Example:

BVH
  Regular rays
     Total intersection tests                     95403
     Interior nodes traversed                   1530971 (16.04 per ray)
     Leaf nodes traversed                        226499 (2.37 per ray)
     Ray/Primitive tests                         337137 (3.53 per ray)
     Intersections found                          94801 (0.99 per ray)
  Shadow rays
     Total intersection tests                      9421
     Interior nodes traversed                    117819 (12.50 per ray)
     Leaf nodes traversed                         15092 (1.60 per ray)
     Ray/Primitive tests                          35331 (3.75 per ray)
     Intersections found                           6381 (0.67 per ray)

*/

#pragma D option quiet

uint64_t kd_intersectionp_tests, kd_intersectionp_hits, kd_intersectionp_primitive_tests;
uint64_t kd_intersectionp_leaf, kd_intersectionp_interior;

:::kdtree_intersectionp_test {
    ++kd_intersectionp_tests;
}

:::kdtree_intersectionp_traversed_leaf_node {
    ++kd_intersectionp_leaf;
}

:::kdtree_intersectionp_traversed_interior_node {
    ++kd_intersectionp_interior;
}

:::kdtree_intersectionp_hit {
    ++kd_intersectionp_hits;
}

:::kdtree_intersectionp_primitive_test {
    ++kd_intersectionp_primitive_tests;
}

uint64_t kd_intersection_tests, kd_intersection_hits, kd_intersection_primitive_tests;
uint64_t kd_intersection_leaf, kd_intersection_interior;

:::kdtree_intersection_test {
    ++kd_intersection_tests;
}

:::kdtree_intersection_traversed_leaf_node {
    ++kd_intersection_leaf;
}

:::kdtree_intersection_traversed_interior_node {
    ++kd_intersection_interior;
}

:::kdtree_intersection_hit {
    ++kd_intersection_hits;
}

:::kdtree_intersection_primitive_test {
    ++kd_intersection_primitive_tests;
}

uint64_t bvh_intersection_tests, bvh_intersection_primitive_tests;
uint64_t bvh_intersection_primitive_hits, bvh_intersection_traversed_interior;
uint64_t bvh_intersection_traversed_leaf;

:::bvh_intersection_started {
    ++bvh_intersection_tests;
}

:::bvh_intersection_traversed_interior_node {
    ++bvh_intersection_traversed_interior;
}

:::bvh_intersection_traversed_leaf_node {
    ++bvh_intersection_traversed_leaf;
}

:::bvh_intersection_primitive_test {
    ++bvh_intersection_primitive_tests;
}

:::bvh_intersection_primitive_hit {
    ++bvh_intersection_primitive_hits;
}

uint64_t bvh_intersectionp_tests, bvh_intersectionp_primitive_tests;
uint64_t bvh_intersectionp_primitive_hits, bvh_intersectionp_traversed_interior;
uint64_t bvh_intersectionp_traversed_leaf;

:::bvh_intersectionp_started {
    ++bvh_intersectionp_tests;
}

:::bvh_intersectionp_traversed_interior_node {
    ++bvh_intersectionp_traversed_interior;
}

:::bvh_intersectionp_traversed_leaf_node {
    ++bvh_intersectionp_traversed_leaf;
}

:::bvh_intersectionp_primitive_test {
    ++bvh_intersectionp_primitive_tests;
}

:::bvh_intersectionp_primitive_hit {
    ++bvh_intersectionp_primitive_hits;
}

/* print it */

dtrace:::END {
    printf("Kd Tree\n");
    printf("  Regular rays\n");
    printf("    Total intersection tests               %12d\n", kd_intersection_tests);
    printf("    Interior nodes traversed               %12d (%d.%02d per ray)\n",
           kd_intersection_interior, RATIO_PARTS(kd_intersection_interior, 
                                                 kd_intersection_tests));
    printf("    Leaf nodes traversed                   %12d (%d.%02d per ray)\n",
           kd_intersection_leaf, RATIO_PARTS(kd_intersection_leaf, 
                                             kd_intersection_tests));
    printf("    Ray/Primitive tests                    %12d (%d.%02d per ray)\n",
           kd_intersection_primitive_tests, 
           RATIO_PARTS(kd_intersection_primitive_tests, 
                       kd_intersection_tests));
    printf("    Intersections found                    %12d (%d.%02d per ray)\n",
           kd_intersection_hits, RATIO_PARTS(kd_intersection_hits,
                                             kd_intersection_tests));
    printf("  Shadow rays\n");
    printf("    Total intersection tests               %12d\n", kd_intersectionp_tests);
    printf("    Interior nodes traversed               %12d (%d.%02d per ray)\n",
           kd_intersection_interior, RATIO_PARTS(kd_intersectionp_interior, 
                                                 kd_intersectionp_tests));
    printf("    Leaf nodes traversed                   %12d (%d.%02d per ray)\n",
           kd_intersectionp_leaf, RATIO_PARTS(kd_intersectionp_leaf, 
                                              kd_intersectionp_tests));
    printf("    Ray/Primitive tests                    %12d (%d.%02d per ray)\n",
           kd_intersectionp_primitive_tests, 
           RATIO_PARTS(kd_intersectionp_primitive_tests, 
                       kd_intersectionp_tests));
    printf("    Intersections found                    %12d (%d.%02d)\n",
           kd_intersectionp_hits, RATIO_PARTS(kd_intersectionp_hits,
                                              kd_intersectionp_tests));
    printf("\n");

    printf("BVH\n");
    printf("  Regular rays\n");
    printf("     Total intersection tests              %12d\n", bvh_intersection_tests);
    printf("     Interior nodes traversed              %12d (%d.%02d per ray)\n",
           bvh_intersection_traversed_interior, 
           RATIO_PARTS(bvh_intersection_traversed_interior, bvh_intersection_tests));
    printf("     Leaf nodes traversed                  %12d (%d.%02d per ray)\n",
           bvh_intersection_traversed_leaf, 
           RATIO_PARTS(bvh_intersection_traversed_leaf, bvh_intersection_tests));
    printf("     Ray/Primitive tests                   %12d (%d.%02d per ray)\n",
           bvh_intersection_primitive_tests,
           RATIO_PARTS(bvh_intersection_primitive_tests, bvh_intersection_tests));
    printf("     Intersections found                   %12d (%d.%02d per ray)\n",
           bvh_intersection_primitive_hits,
           RATIO_PARTS(bvh_intersection_primitive_hits, bvh_intersection_tests));
    printf("  Shadow rays\n");
    printf("     Total intersection tests              %12d\n", bvh_intersectionp_tests);
    printf("     Interior nodes traversed              %12d (%d.%02d per ray)\n",
           bvh_intersectionp_traversed_interior, 
           RATIO_PARTS(bvh_intersectionp_traversed_interior, bvh_intersectionp_tests));
    printf("     Leaf nodes traversed                  %12d (%d.%02d per ray)\n",
           bvh_intersectionp_traversed_leaf, 
           RATIO_PARTS(bvh_intersectionp_traversed_leaf, bvh_intersectionp_tests));
    printf("     Ray/Primitive tests                   %12d (%d.%02d per ray)\n",
           bvh_intersectionp_primitive_tests,
           RATIO_PARTS(bvh_intersectionp_primitive_tests, bvh_intersectionp_tests));
    printf("     Intersections found                   %12d (%d.%02d per ray)\n",
           bvh_intersectionp_primitive_hits,
           RATIO_PARTS(bvh_intersectionp_primitive_hits, bvh_intersectionp_tests));
    printf("\n");
}
