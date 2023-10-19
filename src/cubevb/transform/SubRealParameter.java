package cubevb.transform;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Input.Validate;

@Description("Single dimensional value of a multi-dimensional parameter")
public class SubRealParameter extends RealParameter {
	final public Input<RealParameter> parameterInput = new Input<>("parameter", "parameter that is the source of this sub-parameter", Validate.REQUIRED);
	final public Input<Integer> indexInput = new Input<>("index", "specifies entry of parameter value", Validate.REQUIRED);

	int index;
	RealParameter parameter;
	
	double value, storedValue;

	public SubRealParameter() {
		// we only want one RealParameter and subset specification as input
		lowerValueInput.setRule(Validate.FORBIDDEN);
		upperValueInput.setRule(Validate.FORBIDDEN);
		valuesInput.setRule(Validate.FORBIDDEN);
	}
	
	@Override
	public void initAndValidate() {
		parameter = parameterInput.get();
		index = indexInput.get();
		if (index < 0 || index >= parameter.getDimension()) {
			throw new IllegalArgumentException("index out of range: should be from 0 to " + parameter.getDimension() + " but is " + index);
		}		
	}
	
	@Override
	public Double getValue() {
		value = parameter.getValue(index);
		return value;
	}
	
	@Override
	public Double getValue(int param) {
		return getValue();
	}
	
	
	@Override
	public double getArrayValue() {
		return getValue();
	}

	@Override
	public double getArrayValue(int i) {
		return getValue();
	}

	@Override
	public void setValue(Double value) {
		parameter.setValue(index, value);
	}
	
	@Override
	public void setValue(int i, Double value) {
		setValue(value);
	}
	

	@Override
	public int getDimension() {
		return 1;
	}

	@Override
	public void setDimension(int dim) {
		if (dim == 1) {
			return;
		}
		throw new RuntimeException("Cannot be implemented");
	}

	
	@Override
	public Double getLower() {
		return parameter.getLower();
	}
	
	@Override
	public void setLower(Double lower) {
		parameter.setLower(lower);
	}
	
	@Override
	public Double getUpper() {
		return parameter.getUpper();
	}
	
	@Override
	public void setUpper(Double upper) {
		parameter.setUpper(upper);
	}

	@Override
	public void setBounds(Double lower, Double upper) {
		parameter.setBounds(lower, upper);
	}
	
    @Override
    public void assignFromWithoutID(StateNode other) {
		throw new RuntimeException("Cannot be implemented");
    }
    
    @Override
    public StateNode getCurrent() {
    	return this;
    }

	
	
	
	
	@Override
	protected void store() {
		storedValue = parameter.getValue(index);
		setSomethingIsDirty(false);
	}
	
	@Override
	public void restore() {
		storedValue = parameter.getValue(index);
		setSomethingIsDirty(false);
	}
	
	@Override
	protected boolean requiresRecalculation() {
		if (parameter.getValue(index) != storedValue) {
			setSomethingIsDirty(true);
			return true;
		}
		return false;
	}
}
