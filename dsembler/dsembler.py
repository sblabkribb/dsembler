from flask import Flask, render_template, jsonify, request
app = Flask(__name__)

import requests
from flask_bootstrap import Bootstrap
from input import FileInputForm, InputForm

from pymongo import MongoClient
client = MongoClient('mongodb', 27017, username="root", password="rootpass")
db = client.haseong


app = Flask(__name__)
Bootstrap(app)

app.config['SECRET_KEY'] = '5791628bb0b13ce0c676dfde280ba245'
app.config['UPLOADED_GENE_DEST'] = '/app/uploads'

@app.route('/')
def home():
    default_data = {
            'gene_seq' : "GCCCAAGGCCGCTTGAGCAAATGCTTATGGCGCAGCTGAACGCTGATCTCTAATACTAAAATCACTGCCGTCGATTGATCATTTGGTTGACTTTTGCCAGATACTGAGGCTGGCTATGGGGAGCTGGCGCAGGTGAAAAAACTGCCGATTTTCCCCATGACCCCATCTGGAATCGCCGCCTGCCTTGCGCTATAGCGGCGACCCTGATTTTCCCCATCTAAAAATAAATAGGGGCCTCGCTTACATGCCGATCAAGTACAAGCCTGAAATCCAGCACTCCGATTTCAAGGACCTGACCAACCTGATCCACTTCCAGAGCATGGAAGGCAAGATCTGGCTTGGCGAACAGCGCATGCTGTTGCTGCAGTTTTCAGCGATGGCCAGCTTTCGCCGGGAAATGGTCAATACCCTGGGCATCGAACGCGCCAAGGGCTTGTTCCTGCGCCATGGTTACCAGTCCGGCCTGAAGGATGCCGAACTGGCCAGGAAG",
            'oligomer_size': 50,
            'overlap_size' : 20,
            'optimal_temp' : 60,
            'cluster_size' : 5,
            'cluster_range' : 5
            }
    form = InputForm(request.form, data=default_data)
    forms = FileInputForm()
    return render_template("index.html", form=form, forms=forms)


@app.route('/list_items', methods=['GET'])
def listing():
   items = list(db.dsembler.find({}, {'_id':False}))
   #print(items)
   return jsonify({'result':'success', 'data':items})


@app.route('/insert_item', methods=['POST'])
def insert_item():
    gene_seq = request.form['gene_seq']
    gene_seq_short = gene_seq[:10]+"..."
    doc = {
            'gene_seq' : gene_seq,
            'gene_seq_short' : gene_seq_short,
            'oligomer_size' : request.form['oligomer_size'],
            'overlap_size' : request.form['overlap_size'],
            'optimal_temp' : request.form['optimal_temp'],
            'temp_range' : request.form['temp_range'],
            'cluster_size' : request.form['cluster_size'],
            'cluster_range' : request.form['cluster_range']
            }
    #print(doc)
    db.dsembler.insert_one(doc)
    return jsonify({'result': 'success', 'msg': 'Item inserted'})

if __name__ == '__main__':
    app.run('0.0.0.0',port=5000,debug=True)

