//  PSO_main.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.metaheuristics.singleObjective.particleSwarmOptimization;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

import jmetal.core.Algorithm;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
//import jmetal.metaheuristics.singleObjective.particleSwarmOptimization.objectiffunction.F01_shifted_sphere;
import jmetal.operators.mutation.Mutation;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.singleObjective.Sphere;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.JMException;

/**
 * Class for configuring and running a single-objective PSO algorithm
 */
public class PSO_main {
  public static Logger      logger_ ;      // Logger object
  public static FileHandler fileHandler_ ; // FileHandler object
	
  /**
   * @param args Command line arguments. The first (optional) argument specifies 
   *             the problem to solve.
   * @throws JMException 
   * @throws IOException 
   * @throws SecurityException 
   * Usage: three options
   *      - jmetal.metaheuristics.mocell.MOCell_main
   *      - jmetal.metaheuristics.mocell.MOCell_main problemName
   *      - jmetal.metaheuristics.mocell.MOCell_main problemName ParetoFrontFile
   */
  public static void main(String [] args) 
  		throws JMException, IOException, ClassNotFoundException {
    Problem   problem   ;  // The problem to solve
    Algorithm algorithm ;  // The algorithm to use
    Mutation  mutation  ;  // "Turbulence" operator
    
    QualityIndicator indicators ; // Object to get quality indicators
        
    HashMap  parameters ; // Operator parameters

    // Logger object and file to store log messages
    logger_      = Configuration.logger_ ;
    fileHandler_ = new FileHandler("PSO_main.log"); 
    logger_.addHandler(fileHandler_) ;

    //problem = new Sphere("Real", 20) ;
    //problem = new Easom("Real") ;
     //problem = new Griewank("Real", 10) ;

    //problem = new Sphere("Real", 20); 

    problem = new F01_shifted_sphere(10,-450,10,"Real");
    
     // problem = new F01_shifted_sphere(10, -450,10,"Real");
      //problem = new F02_shifted_schwefel(10, -450,10,"Real");
     // problem = new F03_shifted_rotated_high_cond_elliptic(10, -450,10,"Real");
      // problem= new F04_shifted_schwefel_noise(10, -450,10,"Real");
      //problem= new F05_schwefel_global_opt_bound (10,-310,10,"Real");
      // problem = new F06_shifted_rosenbrock (10,390,10,"Real");
       //problem = new F07_shifted_rotated_griewank (10,-180,10,"Real");
       // problem = new F08_shifted_rotated_ackley_global_opt_bound(10, -140,10,"Real");
     //problem = new F09_shifted_rastrigin (10,-330,10,"Real");
        //problem= new F10_shifted_rotated_rastrigin (10,-330,10,"Real");
   //problem= new F11_shifted_rotated_weierstrass (10,90,10,"Real");
      //problem = new F12_schwefel (10, -460,10,"Real");
    // problem = new F13_shifted_expanded_griewank_rosenbrock (30,-130,10,"Real");
      //problem= new F14_shifted_rotated_expanded_scaffer(10,-300,10,"Real");
      //problem = new F15_hybrid_composition_1 (10,120,10,"Real");
    //problem = new F16_rotated_hybrid_composition_1 (10,120,10,"Real");
       //problem = new F17_rotated_hybrid_composition_1_noise (10,120,10,"Real");
      //problem =   new F18_rotated_hybrid_composition_2 (10,10,10,"Real");
  //problem = new F19_rotated_hybrid_composition_2_narrow_basin_global_opt(10,10,10,"Real");
      //problem = new F20_rotated_hybrid_composition_2_global_opt_bound (10,10,10,"Real");
        //problem = new F21_rotated_hybrid_composition_3 (10,360,10,"Real");
        //problem= new F22_rotated_hybrid_composition_3_high_cond_num_matrix (10,360,10,"Real");
       //problem = new F23_noncontinuous_rotated_hybrid_composition_3 (10,360,10,"Real");
 //problem= new F24_rotated_hybrid_composition_4 (10,260,10,"Real");
          //problem = new F25_rotated_hybrid_composition_4_bound (10,260,10,"Real");

    
    
    
    algorithm = new PSOModify(problem) ;
    
    // Algorithm parameters
    algorithm.setInputParameter("swarmSize",80);
    algorithm.setInputParameter("maxIterations",250);
    int threads = 3; 
    //hh: Passer le nombre de threads en paramètre
    algorithm.setInputParameter("threadNumber",threads);
    int migrationProbability = 10; 
    algorithm.setInputParameter("migrationProbability",migrationProbability);
    //algorithm.setInputParameter("pathIndicators","C:\\Users\\Houda Aba\\Documents\\NetBeansProjects\\jmetall\\src\\jmetal\\problems\\ZDT\\ZDT1.pf");
    algorithm.setInputParameter("pathIndicators","ZDT1.pf");
    
    algorithm.setInputParameter("topology",PSOModify.RING);
    //algorithm.setInputParameter("topology",PSOModify.TORUS);
      //algorithm.setInputParameter("topology",PSOModify.CHAIN);
    parameters = new HashMap() ;
    parameters.put("probability", 1.0/problem.getNumberOfVariables()) ;
    parameters.put("distributionIndex", 20.0) ;
    mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);                    

    algorithm.addOperator("mutation", mutation);

    // Execute the Algorithm 
    PSOModify.initTime = System.currentTimeMillis();
    SolutionSet population = algorithm.execute();
   /* long estimatedTime = System.currentTimeMillis() - initTime;
    
    // Result messages 
    logger_.info("Total execution time: "+estimatedTime + "ms");
    logger_.info("Objectives values have been writen to file FUN");
    population.printObjectivesToFile("FUN");
    logger_.info("Variables values have been writen to file VAR");
    population.printVariablesToFile("VAR");  */                       
  } //main
} // PSO_main
