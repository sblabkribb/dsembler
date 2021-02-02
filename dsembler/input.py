from flask_wtf import FlaskForm
from wtforms import Form, StringField, IntegerField, FloatField, FileField, SubmitField, TextAreaField
from wtforms.validators import DataRequired, NumberRange, InputRequired, EqualTo

class FileInputForm(FlaskForm):
    genes = FileField()

class InputForm(FlaskForm):
    gene_seq = TextAreaField('DNA sequence', render_kw={'rows':5})
    oligomer_size = IntegerField(
            'Oligomer length', 
            validators=[InputRequired(), NumberRange(min=20)],
            render_kw = {'size':10}
            ) #TODO overlap and oligomer must have at least 30 bp difference and oligomer is not less than the overlap length
    overlap_size = IntegerField(
            'Overlap length', 
            validators=[InputRequired(), NumberRange(min=20)]
            )
    optimal_temp = FloatField('Melting Temperature of overlap(ºC):', validators=[InputRequired()])
    temp_range = FloatField('± (ºC):', default=2.5)
    cluster_size = IntegerField('Oligomers per cluster', validators=[InputRequired()])
    cluster_range = IntegerField('±:', validators=[InputRequired(), NumberRange(min=1, max=10)])
