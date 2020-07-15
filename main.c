/**
 * @file cgp.h
 *
 * @brief Cartesian Genetic Programming Algorithm
 *
 * @details This file implements the Cartesian Genetic Programming Algorithm based in 
 * those materials listed below.
 *  -Cartesian Genetic Programming - ISBN-10: 3642269982 ISBN-13: 978-3642269981
 *  -How to evolve complex circuits from scratch - DOI: 10.1109/ICES.2014.7008732
 *  -CGP with Guided and Single Active Mutations for Designing CLCs - DOI: 
 * 
 * @author Lucas Augusto Müller de Souza (lucasmuller@ice.ufjf.br)
 * Computational Engineering student at Universidade Federal de Juiz de Fora
 *
 *
 * @copyright Distributed under the Mozilla Public License 2.0 ( https://opensource.org/licenses/MPL-2.0 )
 *
 * @code available at https://github.com/ciml/ciml-lib/tree/applied-soft-computing-2019
 * @see https://github.com/lucasmullers/
 *
 * Created on: january 15, 2019
 * Updated on: october 27, 2019
 *
 */


#include <time.h>
#include <bdd.h>
#include "include/cgp.h"

int mutation;

/**
* @brief Function that implements the cartesian genetic programming evolutionary process to
* evolve CLCs. The process start with a random population and ends with the first factible
* circuit found.
* @param population - the population struct that will be used into the evolution
* @param table - the table struct that stores the circuits information
* @param gates - the logical gates used into the evolution
* @return 1 if a factible ciruit is found and 0 otherwise
*/
int evolves_cgp_bdd(Individual *population, int *gates)
{
    evaluate_parent_sat_count(population);
    evaluate_population_sat_count(population);

    int best_individual = find_best_individual_sat_count(population);
    finds_individual_active_genes(&population[best_individual]);
    get_max_depth(&population[best_individual]);
    get_num_gates(&population[best_individual]);
    count_num_transistors_individual(&population[best_individual]);
    set_parent(population, best_individual);

    clone_parent(population);
    fprintf(out_file, "--------------------------\n");
    fflush(out_file);
    fprintf(out_file, "Eval.\tIndv.\tGene\tScore\tGates\tTrans.\t");
    for(int i = 0; i < table->num_outputs; i++){
        fprintf(out_file, "DO%d\t", i);
    }
    fprintf(out_file, "Gate\tIn.1\tIn.2\tDepth\t");
    for(int i = 0; i < table->num_inputs; i++){
        fprintf(out_file, "DI%d\t", i);
    }
    fprintf(out_file, "\n");
    while (1)
    {
        if(mutation == 1)
            apply_SAM(population, gates);
        else if(mutation == 2)
            apply_SAM_plus_GAM(population, gates);
        else if(mutation == 3)
            apply_PM(population, gates);
        // evaluate_population_sat_count(population);
    
        
        best_individual = find_best_individual_sat_count(population);
        set_parent(population, best_individual);
		//printf("\n Pop:");

        if (population[0].score == 0)
        {
            fprintf(out_file, "SAT COUNT: %ld INDIVIDUAL: %d EVALUATIONS: %ld\n", population[0].score, best_individual, maxeval);
            fflush(out_file);
            break;
        }
        // if (maxeval % 50000 == 0)
        // {
        //     fprintf(out_file, "SAT COUNT: %ld INDIVIDUAL: %d EVALUATIONS: %ld\n", population[0].score, best_individual, maxeval);
        //     fflush(out_file);
        // }
        if(bdd_getnodenum() >= (int) (0.75 * bdd_getallocnum()))
        {
            bdd_gbc();
        }
        if(maxeval - (NPOP - 1) < 0)
        {
            fprintf(out_file, "SAT COUNT: %ld INDIVIDUAL: %d EVALUATIONS: %ld\n", population[0].score, best_individual, maxeval);
            fflush(out_file);
            return 0;
        }

        clone_parent(population);
    }
    fprintf(out_file, "--------------------------\n");
    fflush(out_file);
    print_post_optimization_data(&population[0]);

    return 1;
}

/**
* @brief Function that implements the cartesian genetic programming evolutionary process to
* optimize CLCs. The process start with a factible population and ends with best factible
* circuit found (the best individual is choosed with respect to the lower number of transistors).
* @param population - the population struct that will be used into the evolution
* @param table - the table struct that stores the circuits information
* @param gates - the logical gates used into the evolution
* @return none
*/
void optimize_circuit(Individual *population, int *gates)
{

    fprintf(out_file,"--------------------------\n");
    fflush(out_file);

    int best_individual = 0;
    // clock_t start = clock();
    while (1)
    {
        if (mutation == 1)
            apply_SAM(population, gates);
        else if (mutation == 2)
            apply_SAM_plus_GAM(population, gates);
        else if (mutation == 3)
            apply_PM(population, gates);
        // evaluate_population_sat_count(population);
		

        clear_population_active_genes(population);
        find_population_active_genes(population);
        best_individual = find_optimized_individual(population);
		/*
		*
		* Change to CGP-RL, update of the occurrence and average matrix.
		*
		*/
		/*
		* End of change
		*/
        set_parent(population, best_individual);

        // if (maxeval % 50000 == 0)
        // {
        //     fprintf(out_file,"NUM TRANSISTORS: %d INDIVIDUAL: %d EVALUATIONS: %ld\n", population[0].num_transistors, best_individual, maxeval);
        //     fflush(out_file);
        // }
        if (bdd_getnodenum() >= (int)(0.75 * bdd_getallocnum()))
        {
            bdd_gbc();
        }
        // if((clock() - start) / (double)CLOCKS_PER_SEC >= 3600.0)
        // {
        //     print_post_optimization_data(&population[0]);
        //     start = clock();
        // }
        if(maxeval - (NPOP - 1) < 0)
        {
            fprintf(out_file,"NUM TRANSISTORS: %d INDIVIDUAL: %d EVALUATIONS: %ld\n", population[0].num_transistors, best_individual, maxeval);
            fflush(out_file);
            break;
        }
        clone_parent(population);
    }
    fprintf(out_file,"--------------------------\n");
    print_post_optimization_data(&population[0]);
}

int main(int argc, char const *argv[])
{
    int semente;
    sscanf(argv[2], "seed=%d", &semente);
    sscanf(argv[3], "ncol=%d", &NCOL);
    sscanf(argv[4], "maxeval=%ld", &maxeval);
    sscanf(argv[5], "mutation=%d", &mutation);
    LB = NCOL/2;
    srand(semente);

    if (argc == 7)
    {
        out_file = fopen(argv[6], "w");
    }
    else if(argc == 8)
    {
        out_file = fopen(argv[7], "w");
    }
    else
    {
        out_file = stdout;
    }

    if (mutation == 1)
        fprintf(out_file, "SAM\n");
    else if (mutation == 2)
        fprintf(out_file, "SAM+GAM\n");
    else if (mutation == 3)
        fprintf(out_file, "PM\n");
    else
    {
        fprintf(out_file, "Mutation value isnt valid!\n");
        exit(1);
    }
    fflush(out_file);
    bdd_init(10000000, 100000);

    Individual *population = (Individual *)malloc(sizeof(Individual) * NPOP);
    int gates[NGATES] = {1, 2, 3, 4, 5, 6, 7, 8};

    table = (Table *)malloc(sizeof(Table));
    table_constructor(argv[1]);

    initialize_population(population, gates);

    clock_t begin, end;
    begin = clock();

    if(argc <= 7)
    {
        if (evolves_cgp_bdd(population, gates))
        {
            optimize_circuit(population, gates);
        }
    }
    else if(argc > 7)
    {
        int aux;
        sscanf(argv[6], "ngates=%d", &aux);
        sow_population(&population[0], argv[1], aux);
        clear_individiual_active_genes(&population[0]);

        calculate_individual_sat_count(&population[0]); 
        if (population[0].score != 0)
        {
            printf("Sow population didn't work!\n");
            exit(1);
        }
        clone_parent(population);
        optimize_circuit(population, gates);
    }

    bdd_done();
    end = clock();
    fprintf(out_file, "TOTAL TIME: %f seconds\n", (end - begin) / (double)CLOCKS_PER_SEC);
    fflush(out_file);

    free(population);
    free(table);

    return 0;
}
