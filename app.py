from flask import Flask, render_template, url_for, flash, redirect, request, send_file, send_from_directory, session, g
from flask_bootstrap import Bootstrap
import sqlite3
from main import Assembly
from input import InputForm
import os.path
import time
import glob

app = Flask(__name__)
Bootstrap(app)

app.config['SECRET_KEY'] = '5791628bb0b13ce0c676dfde280ba245'

conn = sqlite3.connect('database.db')
print('Opened database successfully')

# conn.execute("DROP TABLE IF EXISTS members")
# print("table deleted")

conn.execute('CREATE TABLE IF NOT EXISTS variables (id INTEGER PRIMARY KEY AUTOINCREMENT, user TEXT, short_gene_seq TEXT, gene_seq TEXT, oligomer_size INTEGER, overlap_size INTEGER, melting_temp DECIMAL(3,2), temp_range DECIMAL(1, 2), seq_orientation TEXT)')
print('Table created successfully')

conn.execute('CREATE TABLE IF NOT EXISTS members (member TEXT PRIMARY KEY, name TEXT)')
print('Table created successfully')
conn.close()

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
db_path = os.path.join(BASE_DIR, "database.db")

app.config['SECRET_KEY'] = '5791628bb0b13ce0c676dfde280ba245'

def get_db():
    db = getattr(g, '_database', None)
    if db is None:
        db = g._database = sqlite3.connect(db_path)
    return db

def query_db(query, args=(), one=False):
    cur = get_db().execute(query, args)
    rv = cur.fetchall()
    cur.close()
    return (rv[0] if rv else None) if one else rv


@app.route('/', methods = ['POST', 'GET'])
def index():

    form = InputForm(request.form)
    
    if request.method == 'POST':
        user_option = request.form.get('membercheckbox')
        seq_orientation = request.form.get('seqorientation')
        print(seq_orientation)
        timestamp = time.time_ns()
        gene_seq = (form.gene_seq.data).strip()

        short_gene_seq = gene_seq[:10]+"..."
        oligomer_size = form.oligomer_size.data
        overlap_size = form.overlap_size.data
        melting_temp = form.optimal_temp.data
        temp_range = form.temp_range.data
    
        if user_option == "y":
            user = request.form['user']
            member = query_db('select * from members where member = ?', [user], one=True)
            if member is None:
                flash('No such user', 'danger')
            else:
                flash('Success', 'success')
        else:
            user = f'Guest_{timestamp}'
        
        if seq_orientation == "c":
            seq_orient = seq_orientation
        elif seq_orientation != "c":
            seq_orient = "l"

        a = Assembly(gene_seq, oligomer_size, overlap_size, melting_temp, temp_range, user, seq_orient)
                
        comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, repeats = a.oligomer_design()

        a.output(comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, repeats)        

        with sqlite3.connect(db_path) as con:
            cur = con.cursor()
            cur.execute("INSERT INTO variables (user, short_gene_seq, gene_seq, oligomer_size, overlap_size, melting_temp, temp_range, seq_orientation) VALUES (?,?,?,?,?,?,?,?)",(user, short_gene_seq, gene_seq, oligomer_size, overlap_size, melting_temp, temp_range, seq_orient) )
            
            con.commit()
    else:
        user = None
    
    return render_template("form.html", form=form, user=user)
        

@app.route('/prev-results/<int:record_id>')
def prev_results(record_id):
    record_id = str(record_id)
    
    con = sqlite3.connect(db_path)
    con.row_factory=sqlite3.Row

    cur = con.cursor()
    cur.execute("select * from variables where id = ?", [record_id])

    data = cur.fetchone()

    user = data["user"]
    gene_seq = data["gene_seq"]
    oligomer_size = data["oligomer_size"]
    overlap_size = data['overlap_size']
    melting_temp = data["melting_temp"]
    temp_range = data["temp_range"]
    seq_orientation = data["seq_orientation"]

    a = Assembly(gene_seq, oligomer_size, overlap_size, melting_temp, temp_range, user, seq_orientation)
    
    comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, repeats = a.oligomer_design()

    clusters = a.output(comp_clusters, cluster_five_two_three, cluster_ovr, score, fault, repeats)
    return render_template("user.html", clusters=clusters)

@app.route('/return-excel-file/<user>')
def excel_file(user):
    for filename in glob.glob('/app/output/*.xlsx'):
        if user in filename:
            return send_file(filename, as_attachment=True)

@app.route('/return-fasta-files/<user>')
def fasta_files(user):
    for filename in glob.glob('/app/output/*.fasta'):
        if user in filename:
            return send_file(filename, as_attachment=True)
    

@app.route('/signup/',methods = ['POST', 'GET'])
def signup():
    if request.method == 'POST':
        try:
            member = request.form['username']
            name = request.form['name']
            with sqlite3.connect("database.db") as con:
                cur = con.cursor()
                cur.execute("INSERT INTO members (member, name) VALUES (?, ?)",(member, name))
                con.commit()
                flash("User successfully added", 'success')
        except:
            con.rollback()
            flash("User exists already", 'danger')
        
        finally:
            return redirect(url_for('index'))
            con.close()
    else:
        return redirect(url_for('index') + '#signUpModal')


@app.route('/login/',methods = ['POST', 'GET'])
def login():
    global member
    if request.method == 'POST':
        member = request.form['username']
        user = query_db('select * from members where member = ?', [member], one=True)
        if user is None:
            flash('No such user', 'danger')
            redirect(url_for('index') + '#loginModal')
        else:
            flash(f"Logged in as {member}", 'success')
            session['username'] = member
            return logged_user(member)
        return redirect(url_for('index') + '#loginModal')


@app.route('/user/<member>', methods=["POST", "GET"])
def logged_user(member):
    form = InputForm(request.form)
    if request.method == "POST":
        con = sqlite3.connect(db_path)
        con.row_factory=sqlite3.Row
        
        cur = con.cursor()
        cur.execute("select * from variables where user = ?", [member])

        rows = cur.fetchall()
    return render_template('form.html', form=form, rows=rows, login=True)


@app.route('/logout')
def logout(): 
    session.pop('username', None)
    return redirect(url_for('index'))


@app.route('/list')
def lists():
    con = sqlite3.connect(db_path)
    con.row_factory=sqlite3.Row
    
    cur = con.cursor()
    cur.execute("select * from variables")

    rows = cur.fetchall();
    return render_template("list.html", rows=rows)


@app.route('/list1')
def lists1():
    con = sqlite3.connect(db_path)
    con.row_factory=sqlite3.Row
    
    cur = con.cursor()
    cur.execute("select * from members")

    rows = cur.fetchall();
    return render_template("list copy.html", rowe=rows)


