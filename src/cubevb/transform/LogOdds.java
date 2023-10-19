package cubevb.transform;

public class LogOdds extends Transform {

	@Override
	public double backward(double z) {
		double r = invlogit(z);
		if (parameter.getValue() != r) {
	        parameter.setValue(r);			
		}
		return r;
	}
	
	@Override
	public double forward(double x) {		
		return logit(x);
	}
	
	@Override
	public double jacobian_det(double x) {
		// TODO: is this correct ???
		return 0; 
	}

	@Override
	public double backward(int dim, double oldValue, double newValue) {
		double r = invlogit(newValue);
		if (parameter.getValue(dim) != r) {
	        parameter.setValue(dim, r);			
		}
		return r;
	}

	@Override
	public double jacobian_det(int dim, double z) {
		// TODO: is this correct ???
		return 0;
	}
}
