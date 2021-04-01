from flask import Flask, render_template, url_for, flash, redirect, request, send_file, send_from_directory, session, g
from flask_bootstrap import Bootstrap
import sqlite3
from main import DnaAssemblyDesigner as dad
import os.path
import time

app = Flask(__name__)
Bootstrap(app)

app.config['SECRET_KEY'] = '5791628bb0b13ce0c676dfde280ba245'

conn = sqlite3.connect('database.db')
print('Opened database successfully')

# conn.execute("DROP TABLE IF EXISTS members")
# print("table deleted")

conn.execute('CREATE TABLE IF NOT EXISTS variables (id INTEGER PRIMARY KEY AUTOINCREMENT, user TEXT, short_gene_seq TEXT, gene_seq TEXT, oligomer_size INTEGER, overlap_size INTEGER, melting_temp DECIMAL(3,2), temp_range DECIMAL(1, 2), cluster_size INTEGER, cluster_range INTEGER)')
print('Table created successfully')

conn.execute('CREATE TABLE IF NOT EXISTS members (member TEXT PRIMARY KEY, name TEXT)')
print('Table created successfully')
conn.close()

dad = dad()

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

@app.route('/', methods=['POST', 'GET'])
def index():
    return render_template('form.html')

@app.route('/addrec', methods = ['POST', 'GET'])
def addrec():
    global user, short_gene_seq, gene_seq, oligomer_size, overlap_size, melting_temp, temp_range, cluster_size, cluster_range
    if request.method == 'POST':
        user_option = request.form.get('membercheckbox')
        timestamp = time.time_ns()
        gene_seq = str(request.form['gene_seq']).strip()
        short_gene_seq = gene_seq[:10]+"..."
        oligomer_size = int(request.form['oligomer_size'])
        overlap_size = int(request.form['overlap_size'])
        melting_temp = float(request.form['melting_temp'])
        temp_range = float(request.form['temp_range'])
        cluster_size = int(request.form['cluster_size'])
        cluster_range = int(request.form['cluster_range'])
        if user_option == "y":
            user = request.form['user']
            member = query_db('select * from members where member = ?', [user], one=True)
            if member is None:
                flash('No such user', 'danger')
            else:
                flash('Success', 'success')
        else:
            user = f'Guest_{timestamp}'
        
        clusters = dad.design_oligomers(gene_seq, oligomer_size, overlap_size, melting_temp, temp_range, cluster_size, cluster_range, user)
        
        with sqlite3.connect(db_path) as con:
            cur = con.cursor()
            cur.execute("INSERT INTO variables (user, short_gene_seq, gene_seq, oligomer_size, overlap_size, melting_temp, temp_range, cluster_size, cluster_range) VALUES (?,?,?,?,?,?,?,?,?)",(user, short_gene_seq, gene_seq, oligomer_size, overlap_size, melting_temp, temp_range, cluster_size, cluster_range) )
            
            con.commit()
            if user_option == "y":
                flash("Record successfully added", "success")
    else:
        clusters = None
    return render_template("form.html", user=user)
        

@app.route('/prev-results/<int:record_id>')
def prev_results(record_id):
    global user
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
    cluster_size = data["cluster_size"]
    cluster_range = data["cluster_range"]

    clusters = dad.design_oligomers(gene_seq, oligomer_size, overlap_size, melting_temp, temp_range, cluster_size, cluster_range, user)
    return render_template("user.html", clusters=clusters)

@app.route('/return-excel-file/')
def excel_file():
    return send_file(f'/app/output/oligomers_{user}.xlsx', as_attachment=True)


@app.route('/return-fasta-files/')
def fasta_files():
    return send_file(f'/app/output/fasta_file_{user}.zip', as_attachment=True)


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
    if request.method == "POST":
        con = sqlite3.connect(db_path)
        con.row_factory=sqlite3.Row
        
        cur = con.cursor()
        cur.execute("select * from variables where user = ?", [member])

        rows = cur.fetchall()
    return render_template("form.html", rows=rows, login=True)


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


