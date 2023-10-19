package cubevb.transform;

public class Log extends Transform {

	@Override
	public double forward(double x) {
		double r = Math.log(x); 
		return r;
	}

	@Override
	public double backward(double z) {
		double r = Math.exp(z);
		if (parameter.getValue() != r) {
	        parameter.setValue(r);			
		}
		return r;
	}

	@Override
	public double jacobian_det(double z) {
		return z;
	}

	@Override
	public double backward(int dim, double oldValue, double newValue) {
		double r = Math.exp(newValue);
		if (parameter.getValue(dim) != r) {
	        parameter.setValue(dim, r);			
		}
		return r;
	}

	@Override
	public double jacobian_det(int dim, double z) {
		return z;
	}

}
