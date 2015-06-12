package stats;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.Scanner;

import errors.StatsCalcException;
import tools.SimDist;

public class SimulationParser {
	
	private final int NUM_TESTS = 5;
	
	private String neut_path = "sim_data" + File.separator + "neutral_simulation.tsv";
	private String sel_path = "sim_data" + File.separator + "selection_simulation.tsv";
	
	public SimulationParser(File sim_dir) {
		
		this.neut_path = sim_dir.getAbsolutePath() + File.separator + "neutral_simulation.tsv";
		this.sel_path = sim_dir.getAbsolutePath() + File.separator + "selection_simulation.tsv";
	}
	
	public SimDist[] getNeutralSimulations() throws StatsCalcException {
		
		return parseSimulatedData(neut_path);
		
	}
	
	public SimDist[] getSelectedSimulations() throws StatsCalcException {
		
		return parseSimulatedData(sel_path);
	}
	
	private SimDist[] parseSimulatedData(String file_path) throws StatsCalcException {
		
		SimDist[] dists = createDistributions();
		
		try {
			Scanner scan = new Scanner(new File(file_path));
			scan.nextLine();//skip header line
			
			while(scan.hasNext()) {
				
				String[] vals = scan.nextLine().split("\\s+");
				
				if(vals.length != NUM_TESTS) {
					String err_type = "NumTestError\tIrregularities in number of simulated data columns";
					throw new StatsCalcException(err_type);
				}
				
				//This function depends a lot on the way the sim file is set up
				//AND the SimDist type values (it was engineered this way)
				for(int i = 0; i < NUM_TESTS; i++) {
					double val = Double.parseDouble(vals[i]);
					dists[i].addSimValue(val);
				}
			}
		
		} catch (FileNotFoundException e) {
			String err_type = "IOError\tCould not find simulated data file";
			throw new StatsCalcException(err_type);
		} catch (NumberFormatException e) {
			String err_type = "ParsingError\tInvalid simulated data values";
			throw new StatsCalcException(err_type);
		}
		
		return dists;
	}
	
	/*
	 * Simulations have specified boundaries (according to simulation definition)
	 * 		[0] = iHS bounds [-6,6]
	 * 		[1] = iHH bounds [-3,5]
	 * 		[2] = Fst bounds [-1,6] 
	 * 		[3] = DDAF bounds [-1,1]
	 * 		[4] = XPEHH bounds [-3,8]
	 */
	private SimDist[] createDistributions() {
		
		SimDist[] dists = new SimDist[NUM_TESTS];
		
		dists[SimDist.IHS_TYPE] = new SimDist(-6, 6);
		dists[SimDist.IHH_TYPE] = new SimDist(-3, 5);
		dists[SimDist.FST_TYPE] = new SimDist(-1, 6);
		dists[SimDist.DDAF_TYPE] = new SimDist(-1, 1);
		dists[SimDist.XPEHH_TYPE] = new SimDist(-3, 8);
		
		return dists;
	}
	
	@SuppressWarnings("unused")
	private void printSimDistArr(SimDist[] dists) {
		
		for(int i = 0; i < NUM_TESTS; i++) {
			List<Double> vals = dists[i].getSimVals();
			for(int j = 0; j < vals.size(); j++)
				System.out.print(vals.get(j) + "  ");
			
			System.out.println();
		}
	}

}
