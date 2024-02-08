package cubevb.initialiser;


import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Distribution;
import beast.base.inference.Logger;
import beast.base.inference.Operator;
import beast.base.inference.OperatorSchedule;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beastlabs.inference.SimulatedAnnealing;

@Description("Finds a starting state by running simulated annealing")
public class MAPInitialiserBySimulatedAnnealing extends BEASTObject implements StateNodeInitialiser {
    final public Input<State> startStateInput = new Input<>("state", "elements of the state space");
    final public Input<Distribution> posteriorInput = new Input<>("posterior", "probability distribution to sample over (e.g. a posterior)",
                    Input.Validate.REQUIRED);
    final public Input<List<Operator>> operatorsInput = new Input<>("operator", "operator for generating proposals in MCMC state space",
                    new ArrayList<>());

    final public Input<Integer> saRepeatCountInput = new Input<>("SARepeats", 
    		"number of repeast of running simulated annealing for initialisation", 25);
    final public Input<Integer> chainLengthInput = new Input<>("chainLength", 
    		"number of MCMC samples per simulated annealing run for initialisation. "
    		+ "In total chainLength x SARepeats MCMC steps will be taken", 1000);
    final public Input<Double> startTemp = new Input<Double>("startTemp","starting temperature (default 0.1)", 0.1);
    final public Input<Double> endTemp = new Input<Double>("endTemp","end temperature. Together with startTemp this " +
			"determines the temperature trajectory (default 0.01)", 0.01);

    
    protected State state;
    protected Distribution posterior;

    @Override
	public void initAndValidate() {
		state = startStateInput.get();
		posterior = posteriorInput.get();
	}

	@Override
	public void initStateNodes()  {
		long start = System.currentTimeMillis();
		try {
		
		MySimulatedAnnealing sa = new MySimulatedAnnealing(posterior, state, new OperatorSchedule(), operatorsInput.get(),
				startTemp.get(), 
				endTemp.get(), 
				chainLengthInput.get());

		String bestXML = "";
		double bestLogP = Double.NEGATIVE_INFINITY;
		
		int saRepeatCount = saRepeatCountInput.get();
		for (int i = 0; i < saRepeatCount; i++) {
				Log.info.print(i+"\t");
				sa.doLoop();
			
			if (sa.getLogLikelihood() > bestLogP) {
				bestXML = state.toXML(0);
				bestLogP = sa.getLogLikelihood();
			}

			// state.storeToFile(i);
			// long end = System.currentTimeMillis();
			// Log.info("log posterior = " + bestLogP + " working for " +  (end-start)/1000 + " seconds");
		}		
		
		long end = System.currentTimeMillis();
		Log.info("log posterior = " + bestLogP + " working for " +  (end-start)/1000.0 + " seconds");

		
		state.fromXML(bestXML);
		state.storeToFile(0);
		for (StateNode sn : state.stateNodeInput.get()) {
			if (sn instanceof TreeInterface) {
				Log.info(sn.getID() + ": " + ((TreeInterface) sn).getRoot().toNewick());
			} else {
				Log.info(sn.getID() + ": " + sn.toString());
			}
		}
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

	class MySimulatedAnnealing extends SimulatedAnnealing {
		public MySimulatedAnnealing(Distribution posterior, State state,
				OperatorSchedule operatorSchedule, List<Operator> operators, double startTemp, double endTemp,
				long chainLength) throws IOException {
			Logger logger = new Logger();
			logger.initByName("log", posterior, "logEvery", 1000000);
			logger.init();
			initByName("distribution", posterior, "state", state, "operatorschedule", operatorSchedule,
					"operator", operators, "startTemp", startTemp, "endTemp", endTemp, "chainLength", chainLength,
					"logger", logger);
		}

		protected void doLoop() throws IOException {
			super.doLoop();
		}
		
		double getLogLikelihood() {
			return super.oldLogLikelihood;
		}
	}
	
	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		stateNodes.addAll(state.stateNodeInput.get());
	}

}
