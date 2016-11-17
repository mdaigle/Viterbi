var fs = require('fs');
const assert = require('assert');

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

var diff;
var time = process.hrtime();

for (var i = 1; i <= 10; i++) {
    if (i == 9) {
        diff = process.hrtime(time);
    }

    console.log("Starting round " + i);
    console.log();
    console.log("Emission Parameters:");
    console.log("State 1:");
    console.log(emissions.state1);
    console.log("State 2:");
    console.log(emissions.state2);
    console.log();
    console.log("Transition Parameters:");
    console.log("From State 1:");
    console.log(transitions.state1);
    console.log("From State 2:");
    console.log(transitions.state2);
    console.log();

    var t_count = {
        "state1": { "state1": 0, "state2": 0 },
        "state2": { "state1": 0, "state2": 0 }
    };
    var n_count = {
        "state1": { "A": 0, "C": 0, "G": 0, "T": 0 },
        "state2": { "A": 0, "C": 0, "G": 0, "T": 0 }
    };
    var s_count = { "state1": 0, "state2": 0 };

    var viterbi_path = getViterbiPath(s_count, t_count, n_count);

    if (i == 10) {
        printPathInfo(viterbi_path, -1);
    } else {
        printPathInfo(viterbi_path, 5);
    }

    var total_state_1_transitions = t_count['state1']['state1'] + t_count['state1']['state2'];

    var total_state_2_transitions = t_count['state2']['state1'] + t_count['state2']['state2'];

    var total_state_1_emissions = n_count['state1']['A'] + n_count['state1']['C'] + n_count['state1']['G'] + n_count['state1']['T'];

    var total_state_2_emissions = n_count['state2']['A'] + n_count['state2']['C'] + n_count['state2']['G'] + n_count['state2']['T'];

    // Set begin state probs
    transitions['begin']['state1'] = s_count['state1'] / (s_count['state1'] + s_count['state2']);

    transitions['begin']['state2'] = s_count['state2'] / (s_count['state1'] + s_count['state2']);

    // Set transition probs
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

    console.log("\n\n");
}

console.log("Node 5.10.1");
console.log("Time: " + (diff[0] + 1e-9 * diff[1]) + " seconds");
console.log("Device: MacBook Pro (Retina, 13-inch, Late 2013)");
console.log("Processor: 2.4 GHz Intel Core i5");
console.log("Memory: 8 GB 1600 MHz DDR3");



function getViterbiPath(state_count, transition_count, nucleotide_count) {
    var model = [
        {
            "state1": Math.log(transitions['begin']['state1']),
            "state2": Math.log(transitions['begin']['state2'])
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
        var curr_state_1 = model[model_index - 1]['state1'];
        var curr_state_2 = model[model_index - 1]['state2'];

        // Get the emission probability for nucleotide for each state
        var emission_state_1 = emissions['state1'][nucleotide];
        var emission_state_2 = emissions['state2'][nucleotide];

        if (i == 0) {
            model[model_index]['state1'] = curr_state_1 + Math.log(emission_state_1);
            model[model_index]['state2'] = curr_state_2 + Math.log(emission_state_2);
        } else {
            // Calculate new state 1 probabilities
            var state_1_to_state_1 = curr_state_1 + Math.log(transitions['state1']['state1']) + Math.log(emission_state_1);
            var state_2_to_state_1 = curr_state_2 + Math.log(transitions['state2']['state1']) + Math.log(emission_state_1);

            // Set the definitive state 1 probability to the maximum
            model[model_index]['state1'] = Math.max(state_1_to_state_1, state_2_to_state_1);

            // Calculate new state 2 probabilities
            var state_1_to_state_2 = curr_state_1 + Math.log(transitions['state1']['state2']) + Math.log(emission_state_2);
            var state_2_to_state_2 = curr_state_2 + Math.log(transitions['state2']['state2']) + Math.log(emission_state_2);

            // Set the definitive state 2 probability to the maximum
            model[model_index]['state2'] = Math.max(state_1_to_state_2, state_2_to_state_2);
        }
    }

    // Traceback to get path
    var path = [];

    for (var i = fasta_contents.length; i > 0; i--) {
        var nucleotide = fasta_contents.charAt(i).toUpperCase();
        var curr_nucleotide = fasta_contents.charAt(i - 1).toUpperCase();

        // Replace non ACTG characters
        if (NUCLEOTIDES.indexOf(nucleotide) < 0) {
            nucleotide = "T";
        }
        // Replace non ACTG characters
        if (NUCLEOTIDES.indexOf(curr_nucleotide) < 0) {
            curr_nucleotide = "T";
        }

        /*if (parseInt(nucleotide) < 6) {
            nucleotide = "A";
        } else {
            nucleotide = "C";
        }*/

        if (i == fasta_contents.length) {
            if (model[i]['state1'] >= model[i]['state2']) {
                console.log("Path Log Probability: " + model[i]['state1']);
                path[i] = 'state1';
                state_count['state1']++;
                nucleotide_count['state1'][curr_nucleotide]++;
            } else {
                console.log("Path Log Probability: " + model[i]['state2']);
                path[i] = 'state2';
                state_count['state2']++;
                nucleotide_count['state1'][curr_nucleotide]++;
            }
            continue;
        }

        // Get probabilities of current states
        var curr_state_1 = model[i]['state1'];
        var curr_state_2 = model[i]['state2'];

        // Get probabilities of next states
        var next_state_1 = model[i + 1]['state1'];
        var next_state_2 = model[i + 1]['state2'];

        // Get the emission probability for nucleotide for next states
        var emission_state_1 = emissions['state1'][nucleotide];
        var emission_state_2 = emissions['state2'][nucleotide];

        if (path[i + 1] == 'state1') {
            // Calculate probabilities for state transitions
            var state_1_to_state_1 = curr_state_1 + Math.log(transitions['state1']['state1']) + Math.log(emission_state_1);
            var state_2_to_state_1 = curr_state_2 + Math.log(transitions['state2']['state1']) + Math.log(emission_state_1);

            if (state_1_to_state_1 == next_state_1) {
                path[i] = 'state1';
                state_count['state1']++;
                transition_count['state1']['state1']++;
                nucleotide_count['state1'][curr_nucleotide]++;
            } else if (state_2_to_state_1 == next_state_1) {
                path[i] = 'state2';
                state_count['state2']++;
                transition_count['state2']['state1']++;
                nucleotide_count['state2'][curr_nucleotide]++;
            } else {
                assert(false, "No transition to state 1 gives the expected probability at index " + i);
            }
            continue;
        } else {
            // Calculate probabilities for state transitions
            var state_1_to_state_2 = curr_state_1 + Math.log(transitions['state1']['state2']) + Math.log(emission_state_2);
            var state_2_to_state_2 = curr_state_2 + Math.log(transitions['state2']['state2']) + Math.log(emission_state_2);

            if (state_1_to_state_2 == next_state_2) {
                path[i] = 'state1';
                state_count['state1']++;
                transition_count['state1']['state2']++;
                nucleotide_count['state1'][curr_nucleotide]++;
            } else if (state_2_to_state_2 == next_state_2) {
                path[i] = 'state2';
                state_count['state2']++;
                transition_count['state2']['state2']++;
                nucleotide_count['state2'][curr_nucleotide]++;
            } else {
                assert(false, "No transition to state 2 gives the expected probability at index " + i + "\n" +
                        state_1_to_state_2 + " " + state_2_to_state_2 + " " + next_state_2);
            }
            continue;
        }
    }

    path.shift();

    /*path.forEach(function(elt, i, arr){
        var nucleotide = fasta_contents.charAt(i).toUpperCase();

        // Replace non ACTG characters
        if (NUCLEOTIDES.indexOf(nucleotide) < 0) {
            nucleotide = "T";
        }
        //if (parseInt(nucleotide) < 6) {
        //    nucleotide = "A";
        //} else {
        //    nucleotide = "C";
        //}

        nucleotide_count[elt][nucleotide]++;
    });*/

    /*var output = "";
    path.forEach(function(elt, i, array){
        if (elt == "state1") {
            output += "F";
        } else {
            output += "L";
        }
    });

    console.log(output);*/

    return path;
}

function printPathInfo(path, k) {
    var hits = [];

    var inState2 = false;
    var hit_num = 0;

    path.forEach(function(elt, i, arr){
        if (elt == 'state2' && !inState2) {
            hits[hit_num] = {"start": i + 1, "length": 1};
            inState2 = true;
        } else if (elt == 'state2' && inState2) {
            hits[hit_num]['length']++;
        } else if (inState2){
            hit_num++;
            inState2 = false;
        }
    });

    console.log("Number of Hits: " + hits.length);
    console.log("Lengths and Locations of First k:");

    if (k < 0) {
        k = hits.length;
    }

    for (var i = 0; i < k; i++) {
        if (hits[i] != undefined) {
            console.log(i+1, "start: " + hits[i]['start'], "end: " + (hits[i]['start'] + hits[i]['length'] - 1), "length: " + hits[i]['length']);
        }
    }
}
