package cubevb.transform;

import beast.base.core.Description;

@Description("Transform that leaves the input unaltered. Useful for random walk parameters with infinite range.")
public class NoTransform extends Transform {

	@Override
	public double forward(double x) {
		return x;
	}

	@Override
	public double backward(double z) {
		if (parameter.getValue() != z) {
	        parameter.setValue(z);			
		}
		return z;
	}

	@Override
	public double jacobian_det(double x) {
		return 0;
	}

	@Override
	public double backward(int dim, double oldValue, double newValue) {
		if (parameter.getValue(dim) != newValue) {
	        parameter.setValue(dim, newValue);			
		}
		return newValue;
	}

	@Override
	public double jacobian_det(int dim, double x) {
		return 0;
	}

}
