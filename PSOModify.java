//  PSO.java
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.singleObjective.particleSwarmOptimization.PSO;
import jmetal.operators.selection.BestSolutionSelection;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.ObjectiveComparator;
import jmetal.util.wrapper.XReal;

/**
 * Class implementing a single-objective PSO algorithm
 */
public class PSOModify extends Algorithm {
	private int maxIterations_;
	private int total_particlesSize_; // le nombre des particlues
	private int total_archiveSize_;// the maximum size for the archive
	private int thread_number;// the number of threads

	private int migrationProbability;

	public SolutionSet total_population;

	private Map<Integer, PSOThread> threads;

	private Map<Integer, Solution> numThread_globalBest_map = null;

	private Map<Integer, SolutionSet> numThread_particules_map = null;

	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	String topology;
	static long initTime; // hh: The initial time for executing the
							// algorithm
	static long endTime; // hh: The end time for executing the
							// algorithm

	public static String RING = "RING";
	public static String TORUS = "TORUS";
	public static String CHAIN = "CHAIN";

	public PSOModify(Problem problem) {
		super(problem);
	}

	/**
	 * Initialize all parameter of the algorithm
	 */
	public void initParams() {
		total_particlesSize_ = ((Integer) getInputParameter("swarmSize")).intValue();
		maxIterations_ = ((Integer) getInputParameter("maxIterations")).intValue();
		thread_number = ((Integer) getInputParameter("threadNumber")).intValue();
		migrationProbability = ((Integer) getInputParameter("migrationProbability")).intValue();
		topology = ((String) getInputParameter("topology")).toString();
		threads = new HashMap<Integer, PSOThread>();

		numThread_globalBest_map = new HashMap<Integer, Solution>();
		numThread_particules_map = new HashMap<Integer, SolutionSet>();

		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		try {
			fileHandler_ = new FileHandler("IMOPSO_main_modif.log");
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		logger_.addHandler(fileHandler_);

		total_population = new SolutionSet(total_particlesSize_);

		int particlesSize = (int) total_particlesSize_ / thread_number;
		int archiveSize = (int) total_archiveSize_ / thread_number;
		for (int i = 0; i < thread_number; i++) {
			try {
				PSOThread thread = new PSOThread(this, i + 1, particlesSize, archiveSize);
				threads.put(i + 1, thread);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	} // initParams

	//////////////////
	class PSOThreadAux {
		private int _id;
		private PSOModify _master;

		private int particlesSize_; // le nombre des particlues
		private int iteration_;// the current number of iteration
		private SolutionSet particles_; // the the particles
		private Solution[] localBest_;// the local best solutions found so far
										// for each particles
		private Solution globalBest_;// the global best solution found
		private double[][] speed_;// the speed_ of each particle
		private Operator polynomialMutation_;// a operator for non uniform
												// mutations
		int evaluations_;
		Comparator comparator_;// Comparator object
		Operator findBestSolution_;

		double r1Max_;
		double r1Min_;
		double r2Max_;
		double r2Min_;
		double C1Max_;
		double C1Min_;
		double C2Max_;
		double C2Min_;
		double WMax_;
		double WMin_;
		double ChVel1_;
		double ChVel2_;

		private SolutionSet trueFront_;
		private double deltaMax_[];
		private double deltaMin_[];
		boolean success_;

		private int pourcentage_migration;

		public PSOThreadAux(PSOModify master, int id, int particlesSize, int archiveSize) {
			_master = master;
			_id = id;
			particlesSize_ = particlesSize;
			particles_ = new SolutionSet(particlesSize_);

			r1Max_ = 1.0;
			r1Min_ = 0.0;
			r2Max_ = 1.0;
			r2Min_ = 0.0;
			C1Max_ = 1.5;
			C1Min_ = 1.5;
			C2Max_ = 1.5;
			C2Min_ = 1.5;
			WMax_ = 0.9;
			WMin_ = 0.1;
			ChVel1_ = 1.0;
			ChVel2_ = 1.0;

			comparator_ = new ObjectiveComparator(0); // Single objective
														// comparator
			HashMap parameters; // Operator parameters

			parameters = new HashMap();
			parameters.put("comparator", comparator_);
			findBestSolution_ = new BestSolutionSelection(parameters);

			evaluations_ = 0;

			polynomialMutation_ = operators_.get("mutation");

			iteration_ = 0;

			success_ = false;
			localBest_ = new Solution[particlesSize_];

			// Create the speed_ vector
			speed_ = new double[particlesSize_][problem_.getNumberOfVariables()];

			deltaMax_ = new double[problem_.getNumberOfVariables()];
			deltaMin_ = new double[problem_.getNumberOfVariables()];
			for (int i = 0; i < problem_.getNumberOfVariables(); i++) {
				deltaMax_[i] = (problem_.getUpperLimit(i) - problem_.getLowerLimit(i)) / 2.0;
				deltaMin_[i] = -deltaMax_[i];
			} // for

			numThread_globalBest_map.put(_id, globalBest_);

			pourcentage_migration = (int) particlesSize_ * migrationProbability / 100;
		}

		// Adaptive inertia
		private double inertiaWeight(int iter, int miter, double wmax, double wmin) {
			// return wmax; // - (((wmax-wmin)*(double)iter)/(double)miter);
			return wmax - (((wmax - wmin) * (double) iter) / (double) miter);
		} // inertiaWeight

		// constriction coefficient (M. Clerc)
		private double constrictionCoefficient(double c1, double c2) {
			double rho = c1 + c2;
			// rho = 1.0 ;
			if (rho <= 4) {
				return 1.0;
			} else {
				return 2 / Math.abs((2 - rho - Math.sqrt(Math.pow(rho, 2.0) - 4.0 * rho)));
			}
		} // constrictionCoefficient

		// velocity bounds
		private double velocityConstriction(double v, double[] deltaMax, double[] deltaMin, int variableIndex,
				int particleIndex) throws IOException {

			return v;
			/*
			 * 
			 * //System.out.println("v: " + v + "\tdmax: " + dmax + "\tdmin: " +
			 * dmin) ; double result;
			 * 
			 * double dmax = deltaMax[variableIndex]; double dmin =
			 * deltaMin[variableIndex];
			 * 
			 * result = v;
			 * 
			 * if (v > dmax) { result = dmax; }
			 * 
			 * if (v < dmin) { result = dmin; }
			 * 
			 * return result;
			 */
		} // velocityConstriction

		/**
		 * Update the speed of each particle
		 * 
		 * @throws JMException
		 */
		private void computeSpeed(int iter, int miter) throws JMException, IOException {
			double r1, r2;
			// double W ;
			double C1, C2;
			double wmax, wmin, deltaMax, deltaMin;
			XReal bestGlobal;

			bestGlobal = new XReal(globalBest_);

			for (int i = 0; i < particlesSize_; i++) {
				XReal particle = new XReal(particles_.get(i));
				XReal bestParticle = new XReal(localBest_[i]);

				// int bestIndividual =
				// (Integer)findBestSolution_.execute(particles_) ;

				C1Max_ = 2.5;
				C1Min_ = 1.5;
				C2Max_ = 2.5;
				C2Min_ = 1.5;

				r1 = PseudoRandom.randDouble(r1Min_, r1Max_);
				r2 = PseudoRandom.randDouble(r2Min_, r2Max_);
				C1 = PseudoRandom.randDouble(C1Min_, C1Max_);
				C2 = PseudoRandom.randDouble(C2Min_, C2Max_);
				// W = PseudoRandom.randDouble(WMin_, WMax_);
				//

				WMax_ = 0.9;
				WMin_ = 0.9;
				ChVel1_ = 1.0;
				ChVel2_ = 1.0;

				C1 = 2.5;
				C2 = 1.5;

				wmax = WMax_;
				wmin = WMin_;
				/*
				 * for (int var = 0; var < particle.size(); var++) { //Computing
				 * the velocity of this particle speed_[i][var] =
				 * velocityConstriction(constrictionCoefficient(C1, C2) *
				 * (inertiaWeight(iter, miter, wmax, wmin) * speed_[i][var] + C1
				 * * r1 * (bestParticle.getValue(var) - particle.getValue(var))
				 * + C2 * r2 * (bestGlobal.getValue(var) -
				 * particle.getValue(var))), deltaMax_, deltaMin_, //[var], var,
				 * i); }
				 */
				C1 = 1.5;
				C2 = 1.5;
				double W = 0.9;
				for (int var = 0; var < particle.size(); var++) {
					// Computing the velocity of this particle
					speed_[i][var] = inertiaWeight(iter, miter, wmax, wmin) * speed_[i][var]
							+ C1 * r1 * (bestParticle.getValue(var) - particle.getValue(var))
							+ C2 * r2 * (bestGlobal.getValue(var) - particle.getValue(var));
				}
			}
		} // computeSpeed

		/**
		 * Update the position of each particle
		 * 
		 * @throws JMException
		 */
		private void computeNewPositions() throws JMException {
			for (int i = 0; i < particlesSize_; i++) {
				// Variable[] particle =
				// particles_.get(i).getDecisionVariables();
				XReal particle = new XReal(particles_.get(i));
				// particle.move(speed_[i]);
				for (int var = 0; var < particle.size(); var++) {
					particle.setValue(var, particle.getValue(var) + speed_[i][var]);

					if (particle.getValue(var) < problem_.getLowerLimit(var)) {
						particle.setValue(var, problem_.getLowerLimit(var));
						speed_[i][var] = speed_[i][var] * ChVel1_; //
					}
					if (particle.getValue(var) > problem_.getUpperLimit(var)) {
						particle.setValue(var, problem_.getUpperLimit(var));
						speed_[i][var] = speed_[i][var] * ChVel2_; //
					}

				}
			}
		} // computeNewPositions

		/**
		 * Apply a mutation operator to some particles in the swarm
		 * 
		 * @throws JMException
		 */
		

		/**
		 * Runs of the SMPSO algorithm.
		 * 
		 * @return a <code>SolutionSet</code> that is a set of non dominated
		 *         solutions as a result of the algorithm execution
		 * @throws JMException
		 */
		public SolutionSet execute() throws JMException, ClassNotFoundException {

			success_ = false;
			globalBest_ = null;
			// ->Step 1 (and 3) Create the initial population and evaluate
			for (int i = 0; i < particlesSize_; i++) {
				Solution particle = new Solution(problem_);
				problem_.evaluate(particle);
				evaluations_++;
				particles_.add(particle);
				if ((globalBest_ == null) || (particle.getObjective(0) < globalBest_.getObjective(0))) {
					globalBest_ = new Solution(particle);
					//
					numThread_globalBest_map.put(_id, globalBest_);
				}

			}

			// -> Step2. Initialize the speed_ of each particle to 0
			for (int i = 0; i < particlesSize_; i++) {
				for (int j = 0; j < problem_.getNumberOfVariables(); j++) {
					speed_[i][j] = 0.0;
				}
			}

			// -> Step 6. Initialize the memory of each particle
			for (int i = 0; i < particles_.size(); i++) {
				Solution particle = new Solution(particles_.get(i));
				localBest_[i] = particle;
			}

			// -> Step 7. Iterations ..
			while (iteration_ < maxIterations_) {
				int bestIndividual = (Integer) findBestSolution_.execute(particles_);
				try {
					// Compute the speed_
					computeSpeed(iteration_, maxIterations_);
				} catch (IOException ex) {
					Logger.getLogger(PSO.class.getName()).log(Level.SEVERE, null, ex);
				}

				// Compute the new positions for the particles_
				computeNewPositions();

				// Mutate the particles_
				// mopsoMutation(iteration_, maxIterations_);

				// Evaluate the new particles_ in new positions
				for (int i = 0; i < particles_.size(); i++) {
					Solution particle = particles_.get(i);
					problem_.evaluate(particle);
					evaluations_++;
				}

				// Actualize the memory of this particle
				for (int i = 0; i < particles_.size(); i++) {
					// int flag = comparator_.compare(particles_.get(i),
					// localBest_[i]);
					// if (flag < 0) { // the new particle is best_ than the
					// older remember
					if ((particles_.get(i).getObjective(0) < localBest_[i].getObjective(0))) {
						Solution particle = new Solution(particles_.get(i));
						localBest_[i] = particle;
					} // if
					if ((particles_.get(i).getObjective(0) < globalBest_.getObjective(0))) {
						Solution particle = new Solution(particles_.get(i));
						globalBest_ = particle;
					} // if

				}
				numThread_particules_map.put(_id, particles_);
				threads.get(_id)._aux.particles_ = particles_;
				numThread_globalBest_map.put(_id, globalBest_);
				if (iteration_ == 100) {
					migration(_id, pourcentage_migration);
				}
				iteration_++;
			}

			// Return a population with the best individual
			SolutionSet resultPopulation = new SolutionSet(1);
			resultPopulation.add(particles_.get((Integer) findBestSolution_.execute(particles_)));

			return resultPopulation;
		} // execute

	}

	class PSOThread extends Thread {

		public PSOThreadAux _aux;
		SolutionSet _population;

		public PSOThread(PSOModify master, int id, int particlesSize, int archiveSize) throws Exception {
			_aux = new PSOThreadAux(master, id, particlesSize, archiveSize);
			numThread_particules_map.put(id, _aux.particles_);
		}

		public PSOThreadAux getPSOThreadAux() {
			return _aux;
		}

		public void run() {
			try {
				long initTimeThread = System.currentTimeMillis();
				_population = _aux.execute();
				long estimatedTimeThread = System.currentTimeMillis() - initTimeThread;
				logger_.info("Estimated Time for thread " + _aux._id + " is: " + estimatedTimeThread);
				logger_.info("Variables values have been writen to file VAR");
				_population.printVariablesToFile("VAR_" + _aux._id);
				

				synchronized (this) {
					for (Solution solution : _population.getSolutionsList_()) {
						total_population.add(solution);
					}
				}

				// Print the results
				// logger_.info("Total execution time: "+estimatedTime + "ms");
				/*
				 * synchronized(this) {
				 * 
				 * logger_.info("Variables values have been writen to file VAR"
				 * ); _population.printVariablesToFile("VAR");
				 * logger_.info("Objectives values have been writen to file FUN"
				 * ); _population.printObjectivesToFile("FUN");
				 * 
				 * }
				 */

				if (total_population.size() == thread_number) {
					endTime = System.currentTimeMillis();
					long estimatedTimeForAllAlgorithm = endTime - initTime;
					logger_.info("Estimated Time for all algorithm is " + estimatedTimeForAllAlgorithm);
					logger_.info("-----------------------------------------------------------------");
					// hud: get the best (the minimum) of the number of
					// solutions
					SolutionSet singleSolutionSet = getTheMinimumObjective(total_population);
					logger_.info("Variables values have been writen to file VAR");
					singleSolutionSet.printVariablesToFile("VAR_" + _aux._id);
					logger_.info("Objectives values have been writen to file FUN");
					singleSolutionSet.printObjectivesToFile("FUN");
					// logger_.info("R2 : " + indicators.getR2(population)) ;
				} // if

			} catch (ClassNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (JMException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}

	}

	@Override
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		initParams();

		for (int i = 0; i < thread_number; i++) {
			threads.get(i + 1).start();
		}
		return null;
	}

	public void migration(int numThread, int migrationPourcentage) {
		if (topology.equals(CHAIN)) {
			migrationChainTopology(numThread, migrationPourcentage);
		} else if (topology.equals(RING)) {
			migrationRingTopology(numThread, migrationPourcentage);
		} else if (topology.equals(TORUS)) {
			// migrationTorusTopology(numThread, migrationPourcentage);
		}

	}

	private Solution migrationSolution(Solution solution) {
		Solution newSolution = new Solution(solution);
		return newSolution;
	}

	private void migrationChainTopology(int numThread, int migrationPourcentage) {
		Integer followingThread = null;
		if (numThread == thread_number) {
			followingThread = 1;
		} else {
			followingThread = numThread + 1;
		}
		Solution actual_solution = numThread_globalBest_map.get(numThread);
		// SolutionSet
		// partiulesFollowing=threads.get(followingThread)._aux.particles_;
		SolutionSet partiulesFollowing = numThread_particules_map.get(followingThread);
		SolutionSet partiulesActual = numThread_particules_map.get(numThread);
		
		BestWorstMigration(numThread, partiulesFollowing, partiulesActual, followingThread);

	}

	private void migrationRingTopology(int numThread, int migrationPourcentage) {
		Integer previousThread = null;
		Integer followingThread = null;
		if (numThread == 1) {
			previousThread = thread_number;
			followingThread = numThread + 1;
		} else if (numThread == thread_number) {
			previousThread = numThread - 1;
			followingThread = 1;
		} else {
			previousThread = numThread - 1;
			followingThread = numThread + 1;
		}

		Solution previous_solution = numThread_globalBest_map.get(previousThread);
		Solution following_solution = numThread_globalBest_map.get(followingThread);
		Solution actual_solution = numThread_globalBest_map.get(numThread);

		synchronized (threads.get(numThread)._aux) {
			BestWorstMigration(previousThread, threads.get(numThread)._aux.particles_, threads.get(previousThread)._aux.particles_, numThread);
		}
		synchronized (threads.get(numThread)._aux) {
			BestWorstMigration(followingThread, threads.get(numThread)._aux.particles_, threads.get(followingThread)._aux.particles_, numThread);
		}
		synchronized (threads.get(previousThread)._aux) {
			BestWorstMigration(numThread, threads.get(previousThread)._aux.particles_, threads.get(numThread)._aux.particles_, previousThread);
		}
		synchronized (threads.get(followingThread)._aux) {
			BestWorstMigration(numThread,threads.get(followingThread)._aux.particles_, threads.get(numThread)._aux.particles_, followingThread);
		}

	}

	private void migrationTorusTopology(int numThread, int migrationPourcentage) {
		Integer previousThread = null;
		Integer followingThread = null;
		if (numThread == 1) {
			previousThread = thread_number;
			followingThread = numThread + 1;
		} else if (numThread == thread_number) {
			previousThread = numThread - 1;
			followingThread = 1;
		} else {
			previousThread = numThread - 1;
			followingThread = numThread + 1;
		}
		Solution actual_solution = numThread_globalBest_map.get(numThread);
		if (thread_number > 3) {
			Integer sizeTorus = thread_number - 3;
			List<Integer> listTorus = new ArrayList<Integer>();
			for (int i = 1; i <= thread_number; i++) {
				if ((i != numThread) && (i != previousThread) && (i != followingThread)) {
					listTorus.add(i);
				}
			}
			int index = (int) Math.random() * (listTorus.size() - 1);
			Integer torusThread = listTorus.get(index);
			Solution torus_solution = numThread_globalBest_map.get(torusThread);
			synchronized (threads.get(numThread)._aux) {
				BestWorstMigration(torusThread,threads.get(numThread)._aux.particles_, threads.get(torusThread)._aux.particles_, numThread);
			}
			synchronized (threads.get(torusThread)._aux) {
				BestWorstMigration(numThread,threads.get(torusThread)._aux.particles_, threads.get(numThread)._aux.particles_, torusThread);
			}
		}

		Solution previous_solution = numThread_globalBest_map.get(previousThread);
		Solution following_solution = numThread_globalBest_map.get(followingThread);

		synchronized (threads.get(numThread)._aux) {
			BestWorstMigration(previousThread,threads.get(numThread)._aux.particles_, threads.get(previousThread)._aux.particles_, numThread);
		}
		synchronized (threads.get(numThread)._aux) {
			BestWorstMigration(followingThread,threads.get(numThread)._aux.particles_, threads.get(followingThread)._aux.particles_, numThread);
		}
		synchronized (threads.get(previousThread)._aux) {
			BestWorstMigration(numThread,threads.get(previousThread)._aux.particles_, threads.get(numThread)._aux.particles_, previousThread);
		}
		synchronized (threads.get(followingThread)._aux) {
			BestWorstMigration(numThread,threads.get(followingThread)._aux.particles_, threads.get(numThread)._aux.particles_, followingThread);
		}

	}
	
/**
 * 
 * @param sourceThread: le thread source: le thread qui va migrer leurs bonnes solutions
 * @param particules: les particules de thread target
 * @param sourceParticules: les particules de thread source
 * @param targetThread: le thread target: le thread qui va accepter les bonne solutions
 */
	private synchronized void BestWorstMigration(Integer sourceThread,SolutionSet particules, SolutionSet sourceParticules, Integer targetThread) {
		System.out.println("Debut migration de la thread " + sourceThread + " vers le thread " +
				targetThread);
		//Solution worstSolution = particules.get(0);
		
		//hh: liste des partciules objet de migration
		System.out.println("la liste des particules de thread " + targetThread + " : " + affiche(particules.getSolutionsList_()));
		
		// hh: Trouver toutes les mauvaises particules � �leminer
		List<Solution> worsts = this.getAllWorst(particules);
		
		//hh: liste des partciules � migrer
		System.out.println("la liste des particules � �leminer : " + affiche(worsts));
				
		// hh: Trouver toutes les indexes correspondant aux les plus mauvais particules
		List<Integer> indexes = new ArrayList<Integer>();
		for(int i = 0; i< worsts.size(); i++) {
			Integer index = this.getWorstIndex(particules.getSolutionsList_(), worsts.get(i));
			indexes.add(index);
		}
		

		//hh: liste des particules de l'ile source
		System.out.println("la liste des particules de l'ile source " + sourceThread + " : "+ affiche(sourceParticules.getSolutionsList_()));
		List<Solution> bests = this.getAllBest(sourceParticules);
		//hh: liste des meilleurs particules de l'ile source
		System.out.println("la liste des meilleurs particules de source � migrer : " + affiche(bests));
		// hh migration
            
               	int pourcentage_migration = (int) sourceParticules.getSolutionsList_().size() * migrationProbability / 100;	
                
		//ArrayList<Double> targetObjectifsList = new ArrayList<Double>();
		//Collections.addAll(targetObjectifsList, );
		for(int i =1; i<=pourcentage_migration; i++) {
			threads.get(targetThread)._aux.particles_.getSolutionsList_().remove(worsts.get(i-1));
			threads.get(targetThread)._aux.particles_.getSolutionsList_().add(indexes.get(i-1), bests.get(i-1));

		}
                
		//hh: liste des particules apr�s migration
		System.out.println("la liste des particules de le thread " + targetThread + " apr�s migration: " + affiche(threads.get(targetThread)._aux.particles_.getSolutionsList_()));
		
	}
	
	private String affiche (List<Solution> solutions) {
		String target = "";
		for(Solution element: solutions) {
			target = target + element.getObjective_()[0] + " - ";
		}
		return target;
	}
	
	private String affiche1 (List<Double> liste) {
		String target = "";
		for(Double element: liste) {
			target = target + element + " - ";
		}
		return target;
	}

	private Solution getWorst(List<Solution> list) {
		//double[] objectives = solution.getObjective_();
		Solution worst = list.get(0);
		for (int i = 0; i < list.size(); i++) {
			if (list.get(i).getObjective_()[0] > worst.getObjective_()[0]) {
				worst = list.get(i);
			}
		}
		return worst;
	}
	
	// hh: Trouver la meilleure solution parmi les particules
	private Solution getBest(ArrayList<Solution> list) {
		Solution best = list.get(0);
		for (int i = 0; i < list.size(); i++) {
			if (list.get(i).getObjective_()[0] < best.getObjective_()[0]) {
				best = list.get(i);
			}
		}
		return best;
	}

	// hh: Trouver les diff�rents particules � �liminer suivant la
	// pourcentage
	private List<Solution> getAllWorst(SolutionSet particules) {
		int pourcentage_migration = (int) particules.getSolutionsList_().size() * migrationProbability / 100;
		List<Solution> list = new ArrayList<Solution>();
		//ArrayList<Solution> objectifsList = (ArrayList<Solution>) particules.getSolutionsList_();
		ArrayList<Solution> objectifsList = new ArrayList<Solution>();
		for(Solution solution : particules.getSolutionsList_()) {
			objectifsList.add(solution);
		}
		
		for (int i = 0; i < pourcentage_migration; i++) {
			Solution worst = this.getWorst(objectifsList);
			
			//Collections.addAll(objectifsList, solution.getObjective_());
			//ArrayList<Double> objectifsList = (ArrayList<Double>) Arrays.asList(solution.getObjective_());
			list.add(worst);
			objectifsList.remove(worst);
			
		}
		
		return list;
	}
	
	// hh: Trouver les meilluers partciules suivant la pourcentage
	private List<Solution> getAllBest(SolutionSet particules) {
		int pourcentage_migration = (int) particules.getSolutionsList_().size() * migrationProbability / 100;
		List<Solution> list = new ArrayList<Solution>();
		
		ArrayList<Solution> objectifsList = new ArrayList<Solution>();
		for(Solution solution : particules.getSolutionsList_()) {
			objectifsList.add(solution);
		}
		
		for (int i = 0; i < pourcentage_migration; i++) {
			Solution best = this.getBest(objectifsList);
			
			//Collections.addAll(objectifsList, solution.getObjective_());
			//ArrayList<Double> objectifsList = (ArrayList<Double>) Arrays.asList(solution.getObjective_());
			list.add(best);
			objectifsList.remove(best);
			
		}
		return list;
	}
	//hh: trouver l'index de l'�lement � �leminer
	private int getWorstIndex(List<Solution> solutions, Solution worst) {
		for (int i=0; i<solutions.size(); i++) {
			Solution solution = solutions.get(i);
			if (Double.compare(solution.getObjective_()[0], worst.getObjective_()[0]) == 0){
				return i;
			}
		}
		return -1;
	}

	private SolutionSet getTheMinimumObjective(SolutionSet total_population) {
		SolutionSet solutionSet = new SolutionSet(1);
		double minimumObjective = total_population.get(0).getObjective(0);
		solutionSet.add(total_population.get(0));
		for (int i = 0; i < total_population.size(); i++) {
			double objective = total_population.get(i).getObjective(0);
			if (objective < minimumObjective) {
				solutionSet.remove(0);
				solutionSet.add(total_population.get(i));
			}
		}
		return solutionSet;
	}

} // PSO
