from flask_wtf import FlaskForm
from wtforms import Form, StringField, IntegerField, FloatField, FileField, SubmitField, TextAreaField, SelectField
from wtforms.validators import DataRequired, NumberRange, InputRequired, EqualTo

class FileInputForm(FlaskForm):
    genes = FileField()

class InputForm(FlaskForm):
    gene_seq = TextAreaField('Gene Sequence:')
    oligomer_size = IntegerField('Target Oligomer Size', validators=[InputRequired(), NumberRange(min=20)]) #TODO overlap and oligomer must have at least 30 bp difference and oligomer is not less than the overlap length
    overlap_size = IntegerField('Target Overlap Size:', validators=[InputRequired(), NumberRange(min=20)])
    optimal_temp = FloatField('Melting Temperature of overlap(ºC):', validators=[InputRequired()])
    temp_range = FloatField('± (ºC):', default=2.5)
    cluster_size = IntegerField('Target Cluster Size:', validators=[InputRequired()])
    cluster_range = IntegerField('±:', validators=[InputRequired(), NumberRange(min=1, max=10)])