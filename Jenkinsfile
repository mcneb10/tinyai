pipeline {
    agent { label 'linux' }
    stages {
        stage('Setup') {
            steps {
                echo 'Setting up...'
                echo '- Updating'
                sh 'sudo apt-get update'
                echo '- Installing build tools'
                sh 'sudo apt-get install -y build-essential perl graphviz'
                sh 'cpan Object::Pad File::Slurper'
            }
        }
        stage('Build') {
            steps {
                echo 'Building...'
                sh 'make -C experiments/xor/'
            }
        }
        stage('Test') {
            steps {
                echo 'Testing...'
                // Do tests
                sh './experiments/xor/xor'
                sh './experiments/xor/xor_independent'
                sh './experiments/xor/xor_changeable_activation_aggregation'
                sh './experiments/xor/xor_changeable_activation_aggregation_independent'

                // Test bonus script
                sh 'perl ./bonus/neuralnet2dot.pl -p ./experiments/xor/fit = 200 ./graph200.gv'
                sh 'perl ./bonus/neuralnet2dot.pl -p ./experiments/xor/fitc = 200 ./graphc200.gv'

                // Test graphviz output
                sh 'dot -Tps ./graph200.gv -o graph200.ps'
                sh 'dot -Tps ./graphc200.gv -o graphc200.ps'
            }
        }
    }
}