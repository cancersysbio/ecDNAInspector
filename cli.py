# cli.py

import argparse
import sys
from ecDNAInspector.cycle_select_pipeline import (
    cycle_data_calculations,
    cycle_clustering,
    visualize_clustering,
    assign_cycle_confidence_by_cluster,
    assign_cycle_confidence_by_selection,
    intrasample_filtering)
from ecDNAInspector.utils import load_config


def main():
    parser = argparse.ArgumentParser(description="Run ecDNAInspector")
    parser.add_argument(
        "-c", "--config", type=str,
        default="configs/default_config.yaml",
        help="Path to YAML configuration file"
    )

    parser.add_argument("--run-metric-calc", action="store_true", help="Run cycle metric calculations")
    parser.add_argument("--run-cluster", action="store_true", help="Run clustering step")
    parser.add_argument("--run-intrasample-filter", action="store_true", help="Run intrasample filtering")
    parser.add_argument("--run-conf-assignment", action="store_true", help="Run confidence assignment")

    parser.add_argument(
        "-conf_type", "--confidence_assignment_type", type=str,
        choices=["by_selection", "by_cluster"],
        help="Type of cycle confidence assignment process"
    )

    parser.add_argument(
        "--high_conf_clusters", nargs='+', type=str, default=[],
        help="List of manually selected high confidence clusters"
    )

    parser.add_argument(
        "--med_conf_clusters", nargs='+', type=str, default=[],
        help="List of manually selected medium confidence clusters"
    )

    parser.add_argument(
        "--low_conf_clusters", nargs='+', type=str, default=[],
        help="List of manually selected low confidence clusters"
    )

    parser.add_argument(
        "--visualize_cluster",
        choices=["unfiltered", "intrasample_filtered"],
        help="Type of cluster visualization to produce"
    )

    args = parser.parse_args()

    config_arg = args.config
    config = load_config(args.config)

    if args.run_metric_calc:
        print("Working on cycle metric calculations...")
        cycle_data_calculations(config_arg)
        print("Cycle metric calculations complete!")

    if args.run_cluster:
        print("Working on clustering...")
        cycle_clustering(config_arg)
        print("Clustering complete!")

    if args.run_conf_assignment:
        print("Working on cycle confidence assignment...")
        if not args.confidence_assignment_type:
            print("Error: --confidence_assignment_type is required when --run-conf-assignment is set.", file=sys.stderr)
            sys.exit(1)

        if args.confidence_assignment_type == "by_cluster":
            if not args.high_conf_clusters and not args.med_conf_clusters and not args.low_conf_clusters:
                print("Error: At least one of --high_conf_clusters/med_conf_clusters/low_conf_clusters are required when --confidence_assignment_type is by_cluster.",
                      file=sys.stderr)
                sys.exit(1)
            else:
                assign_cycle_confidence_by_cluster(config_arg, args.high_conf_clusters, args.med_conf_clusters, args.low_conf_clusters)
                print("Cycle confidence assignment complete!")

        elif args.confidence_assignment_type == "by_selection":
            assign_cycle_confidence_by_selection(config_arg)
            print("Cycle confidence assignment complete!")

        else:
            print(
                "Error: invalid --confidence_assignment_type.",
                file=sys.stderr)
            sys.exit(1)

    if args.run_intrasample_filter:
        print("Please confirm that you have already run the confidence assignment!\nThis is required before intrasample filtering.")
        response = input("Press ENTER to continue, or type 'n' to cancel: ").strip().lower()
        if response == 'n':
            print("Aborting intrasample filtering.")
            sys.exit(0)
        print("Working on intrasample filtering...")
        intrasample_filtering(config_arg)
        print("Intrasample filtering complete!")

    if args.visualize_cluster == "unfiltered":
        print("Working on clustering visualization...")
        cluster_table = config['data']['cycle_level_data_clustered_table']
        cluster_image = config['files']['unfiltered_cluster_image_output']
        visualize_clustering(cluster_table, cluster_image)
        print("Clustering visualization complete!")

    if args.visualize_cluster == "intrasample_filtered":
        print("Working on clustering visualization...")
        cluster_table = config['data']['intrasample_filt_cycle_level_data_w_conf_table']
        cluster_image = config['files']['intrasamp_filtered_cluster_image_output']
        visualize_clustering(cluster_table, cluster_image)
        print("Clustering visualization complete!")


if __name__ == "__main__":
    main()
