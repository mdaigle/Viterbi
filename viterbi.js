var fs = require('fs');
const assert = require('assert');
var math = require('mathjs');

math.config({
    number: 'BigNumber',
    precision: 10
});

var STOP_CODONS = ["TAA", "TAG", "TGA"];
var NUCLEOTIDES = ["A", "C", "G", "T"];

var emissions = {
    "state1": { "A": .25, "C": .25, "G": .25, "T": .25 },
    "state2": { "A": .20, "C": .30, "G": .30, "T": .20 }
}

/*var emissions = {
    "state1": { "A": 1/6, "C": 1/6},
    "state2": { "A": 1/10, "C": 1/2}
}*/

var transitions = {
    "begin": { "state1": .9999, "state2": .0001 },
    "state1": { "state1": .9999, "state2": .0001 },
    "state2": { "state1": .01, "state2": .99}
}

/*var transitions = {
    "begin": { "state1": .95, "state2": .05 },
    "state1": { "state1": .95, "state2": .05 },
    "state2": { "state1": .1, "state2": .9}
}*/

var args = process.argv.slice(2);
var fasta_file = args[0];
//var gbb_file = args[1];
//var summary = args[2];
var fasta_contents = fs.readFileSync(fasta_file, 'utf8');
//var gbb_contents = fs.readFileSync(gbb_file, 'utf8');

// Get interesting stuff from fasta file
fasta_contents = fasta_contents.slice(fasta_contents.indexOf("\n") + 1); // Remove comment line
fasta_contents = fasta_contents.split(">", 1)[0]; // Remove plasmids
fasta_contents = fasta_contents.replace(/\s/g, ""); // Remove whitespace

// Parse Genbank file
/*gbb_contents = gbb_contents.split("ORIGIN")[0];
gbb_contents = gbb_contents.split(/\n/g); // Break into lines
var genbank_sequences = [];
gbb_contents.forEach(function(line){
    var regex = /\s+CDS\s+([0-9]+)\.\.([0-9]+)/g;
    var sequence = regex.exec(line);
    if (sequence && sequence.length == 3) {
        genbank_sequences.push({
            "start": sequence[1],
            "end": sequence[2]
        });
    }
});*/

for (var i = 0; i < 10; i++) {
    var t_count = {
        "state1": { "state1": 0, "state2": 0 },
        "state2": { "state1": 0, "state2": 0 }
    }
    var n_count = {
        "state1": { "A": 0, "C": 0, "G": 0, "T": 0 },
        "state2": { "A": 0, "C": 0, "G": 0, "T": 0 }
    }
    var s_count = { "state1": 0, "state2": 0 }

    var viterbi_path = getViterbiPath(s_count, t_count, n_count);

    console.log(s_count);
    console.log(t_count);
    console.log(n_count);

    var total_state_1_transitions = t_count['state1']['state1'] + t_count['state1']['state2'];

    var total_state_2_transitions = t_count['state2']['state1'] + t_count['state2']['state2'];

    var total_state_1_emissions = n_count['state1']['A'] + n_count['state1']['C'] + n_count['state1']['G'] + n_count['state1']['T'];

    var total_state_2_emissions = n_count['state2']['A'] + n_count['state2']['C'] + n_count['state2']['G'] + n_count['state2']['T'];

    // Set begin state probs
    transitions['begin']['state1'] = s_count['state1'] / (s_count['state1'] + s_count['state2']);

    transitions['begin']['state2'] = s_count['state2'] / (s_count['state1'] + s_count['state2']);

    // Set transision probs
    transitions['state1']['state1'] = t_count['state1']['state1'] / total_state_1_transitions;
    transitions['state1']['state2'] = t_count['state1']['state2'] / total_state_1_transitions;
    transitions['state2']['state1'] = t_count['state2']['state1'] / total_state_2_transitions;
    transitions['state2']['state2'] = t_count['state2']['state2'] / total_state_2_transitions;

    // Set emission probs
    emissions['state1']['A'] = n_count['state1']['A'] / total_state_1_emissions;
    emissions['state1']['C'] = n_count['state1']['C'] / total_state_1_emissions;
    emissions['state1']['G'] = n_count['state1']['G'] / total_state_1_emissions;
    emissions['state1']['T'] = n_count['state1']['T'] / total_state_1_emissions;

    emissions['state2']['A'] = n_count['state2']['A'] / total_state_2_emissions;
    emissions['state2']['C'] = n_count['state2']['C'] / total_state_2_emissions;
    emissions['state2']['G'] = n_count['state2']['G'] / total_state_2_emissions;
    emissions['state2']['T'] = n_count['state2']['T'] / total_state_2_emissions;

    console.log(transitions);
    console.log(emissions);
}

// TODO: use log probabilities

function getViterbiPath(state_count, transition_count, nucleotide_count) {
    var model = [
        {
            "state1": math.log(math.bignumber(transitions['begin']['state1'])),
            "state2": math.log(math.bignumber(transitions['begin']['state2']))
        }
    ];

    for (var i = 0; i < fasta_contents.length; i++) {
        var nucleotide = fasta_contents.charAt(i).toUpperCase();

        // Replace non ACTG characters
        if (NUCLEOTIDES.indexOf(nucleotide) < 0) {
            nucleotide = "T";
        }

        /*if (parseInt(nucleotide) < 6) {
            nucleotide = "A";
        } else {
            nucleotide = "C";
        }*/

        var model_index = i + 1;
        model[model_index] = {};

        // Get probabilities of previous states
        var prev_state_1 = model[model_index - 1]['state1'];
        var prev_state_2 = model[model_index - 1]['state2'];

        // Get the emission probability for nucleotide for each state
        var emission_state_1 = emissions['state1'][nucleotide];
        var emission_state_2 = emissions['state2'][nucleotide];

        // Calculate new state 1 probabilities
        var state_1_to_state_1 = math.add(math.add(prev_state_1, math.log(math.bignumber(transitions['state1']['state1']))), math.log(math.bignumber(emission_state_1)));
        var state_2_to_state_1 = math.add(math.add(prev_state_2, math.log(math.bignumber(transitions['state2']['state1']))), math.log(math.bignumber(emission_state_1)));

        // Set the definitive state 1 probability to the maximum
        model[model_index]['state1'] = math.max(state_1_to_state_1, state_2_to_state_1);

        // Calculate new state 2 probabilities
        var state_1_to_state_2 = math.add(math.add(prev_state_1, math.log(math.bignumber(transitions['state1']['state2']))), math.log(math.bignumber(emission_state_2)));
        var state_2_to_state_2 = math.add(math.add(prev_state_2, math.log(math.bignumber(transitions['state2']['state2']))), math.log(math.bignumber(emission_state_2)));

        // Set the definitive state 2 probability to the maximum
        model[model_index]['state2'] = math.max(state_1_to_state_2, state_2_to_state_2);
    }

    //console.log(model);

    // Traceback to get path
    var path = [];

    for (var i = fasta_contents.length - 1; i >= 0; i--) {
        var model_index = i + 1;

        if (model_index == fasta_contents.length) {
            if (math.largerEq(model[model_index]['state1'], model[model_index]['state2'])) {
                path[i] = 'state1';
                state_count['state1']++;
            } else {
                path[i] = 'state2';
                state_count['state2']++;
            }
            continue;
        }

        var nucleotide = fasta_contents.charAt(i).toUpperCase();

        // Replace non ACTG characters
        if (NUCLEOTIDES.indexOf(nucleotide) < 0) {
            nucleotide = "T";
        }
        /*if (parseInt(nucleotide) < 6) {
            nucleotide = "A";
        } else {
            nucleotide = "C";
        }*/

        // Get probabilities of previous states
        var prev_state_1 = model[model_index - 1]['state1'];
        var prev_state_2 = model[model_index - 1]['state2'];

        // Get probabilities of current states
        var curr_state_1 = model[model_index]['state1'];
        var curr_state_2 = model[model_index]['state2'];

        // Get the emission probability for nucleotide for each state
        var emission_state_1 = emissions['state1'][nucleotide];
        var emission_state_2 = emissions['state2'][nucleotide];

        if (path[i + 1] == 'state1') {
            // Calculate probabilities for state transitions
            var state_1_to_state_1 = prev_state_1 * transitions['state1']['state1'] * emission_state_1;
            var state_2_to_state_1 = prev_state_2 * transitions['state2']['state1'] * emission_state_1;

            if (state_1_to_state_1 == curr_state_1) {
                path[i] = 'state1';
                state_count['state1']++;
                transition_count['state1']['state1']++;
            } else if (state_2_to_state_1 == curr_state_1) {
                path[i] = 'state2';
                state_count['state2']++;
                transition_count['state2']['state1']++;
            } else {
                assert(false, "No transition gives the expected probability");
            }
            continue;
        } else {
            // Calculate probabilities for state transitions
            var state_1_to_state_2 = prev_state_1 * transitions['state1']['state2'] * emission_state_2;
            var state_2_to_state_2 = prev_state_2 * transitions['state2']['state2'] * emission_state_2;

            if (state_1_to_state_2 == curr_state_2) {
                path[i] = 'state1';
                state_count['state1']++;
                transition_count['state1']['state2']++;
            } else if (state_2_to_state_2 == curr_state_2) {
                path[i] = 'state2';
                state_count['state2']++;
                transition_count['state2']['state2']++;
            } else {
                assert(false, "No transition gives the expected probability");
            }
            continue;
        }
    }

    // TODO: Maybe do this inside the traceback loop?
    path.forEach(function(elt, i, arr){
        var nucleotide = fasta_contents.charAt(i).toUpperCase();

        // Replace non ACTG characters
        if (NUCLEOTIDES.indexOf(nucleotide) < 0) {
            nucleotide = "T";
        }
        /*if (parseInt(nucleotide) < 6) {
            nucleotide = "A";
        } else {
            nucleotide = "C";
        }*/

        nucleotide_count[elt][nucleotide]++;
    });

    var output = "";
    path.forEach(function(elt, i, array){
        if (elt == "state1") {
            output += "F";
        } else {
            output += "L";
        }
    });

    console.log(output);

    return path;
}

// Print out path stats
/*var output = "";
path.forEach(function(elt, i, array){
    if (elt == "state1") {
        output += "F";
    } else {
        output += "L";
    }
});

console.log(output);*/

// Print out
